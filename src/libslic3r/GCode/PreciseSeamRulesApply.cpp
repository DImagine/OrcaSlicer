// Precise Seam layer rules applied to a print. Design: docs/HLSD/precise-seam.md
#include "PreciseSeamRulesApply.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include "libslic3r/I18N.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "PreciseSeam.hpp"
#include "PreciseSeamRules.hpp"

namespace Slic3r {
namespace PreciseSeam {

bool uses_precise_seam_rules(const ModelObject &object)
{
    return std::any_of(object.volumes.begin(), object.volumes.end(), [](const ModelVolume *volume) {
        return volume->is_precise_seam() && (!volume->precise_seam_rules.empty() || has_rule_marker(volume->name));
    });
}

void init_precise_seam_data_by_rules(
    std::vector<const ModelVolume*> &strong_volumes_out,
    std::vector<const ModelVolume*> &weak_volumes_out,
    bool &has_strong_out,
    const ModelObject &model_object,
    const LayerNumbering &numbering,
    const ConstLayerPtrsAdaptor &model_layers,
    EffectiveModifiers &effective_out,
    std::vector<ModifierDiagnostic> &diagnostics_out)
{
    strong_volumes_out.clear();
    weak_volumes_out.clear();
    effective_out.clear();
    has_strong_out = false;

    // Rules are parsed once per volume, before the parallel layer processing.
    std::vector<const ModelVolume*> volumes;
    std::vector<LayerRules>         rules;
    for (const ModelVolume* volume : model_object.volumes)
        if (volume->is_precise_seam()) {
            volumes.push_back(volume);
            rules.push_back(parse_layer_rules(volume->name, volume->precise_seam_rules));
        }
    if (volumes.empty())
        return;
    std::vector<RuledModifier> ruled;
    for (size_t i = 0; i < volumes.size(); ++i)
        ruled.push_back({volumes[i]->type(), &rules[i]});
    const EffectiveOrder order = effective_order(ruled);

    // Diagnostics of one volume in source order (the name as a whole, its blocks, then the field) and
    // element order; parsing and empty-set diagnostics of one element keep their order.
    const auto source_rank = [](const RuleDiagnostic &d) {
        return std::pair<size_t, size_t>(d.source.field ? std::numeric_limits<size_t>::max() : d.source.name_block, d.element);
    };
    const auto contains = [](const std::vector<size_t> &list, size_t i) { return std::find(list.begin(), list.end(), i) != list.end(); };

    // Layers of each volume, with its diagnostics in volume order. An empty row means every layer;
    // nullopt leaves the volume out: switched off, or its rules exclude every layer of this object.
    const std::vector<size_t> numbered = numbering.numbered_layers();
    // Height bands of the numbered model layers, matched by the original layer index (numbered_layers()
    // drops that index); built only when some volume selects by height.
    std::vector<LayerBand> bands;
    const bool uses_heights = std::any_of(rules.begin(), rules.end(), [](const LayerRules &r) {
        return std::any_of(r.operations.begin(), r.operations.end(),
                           [](const RuleLayerOperation &op) { return op.set.kind == RuleLayerSet::Kind::Height; });
    });
    if (uses_heights) {
        // Both come from the same print object, so they always have one entry per layer.
        assert(numbering.object_layers.size() == model_layers.size());
        for (size_t layer = 0; layer < numbering.object_layers.size(); ++layer)
            if (numbering.object_layers[layer] != 0)
                bands.push_back({numbering.object_layers[layer], model_layers[layer]->bottom_z(), model_layers[layer]->print_z});
    }
    std::vector<std::optional<std::vector<bool>>> layers(volumes.size());
    for (size_t i = 0; i < volumes.size(); ++i) {
        // Nothing reads the diagnostics of the rules after this point.
        std::vector<RuleDiagnostic> own = std::move(rules[i].diagnostics);
        if (!rules[i].off && !rules[i].restricts_layers())
            layers[i].emplace();
        else if (!rules[i].off) {
            LayerSelection selection = select_layers(rules[i], numbering.count, numbered, bands);
            own.insert(own.end(), std::make_move_iterator(selection.diagnostics.begin()), std::make_move_iterator(selection.diagnostics.end()));
            std::vector<bool> row(numbering.object_layers.size(), false);
            bool              any = false;
            for (size_t layer = 0; layer < row.size(); ++layer) {
                row[layer] = selection.contains(numbering.object_layers[layer]);
                any        = any || row[layer];
            }
            if (any)
                layers[i] = std::move(row);
        }
        std::stable_sort(own.begin(), own.end(),
                         [&source_rank](const RuleDiagnostic &a, const RuleDiagnostic &b) { return source_rank(a) < source_rank(b); });
        // Differences from the object list concern the volume as a whole and follow its own diagnostics.
        if (contains(order.type_changed, i))
            own.push_back({RuleDiagnosticCode::TypeDiffersFromTree, {}, 0, {}});
        if (contains(order.reordered, i))
            own.push_back({RuleDiagnosticCode::OrderDiffersFromTree, {}, 0, {}});
        for (RuleDiagnostic &diagnostic : own)
            diagnostics_out.push_back({volumes[i], std::move(diagnostic)});
        if (rules[i].name_overridden)
            BOOST_LOG_TRIVIAL(info) << "[PreciseSeamRules] override: object_id=" << model_object.id().id
                                    << " volume_id=" << volumes[i]->id().id << " the rules field replaces the name rules";
        // nowarn hides the differences from the user, not from the log.
        if (contains(order.suppressed, i))
            // IDs only: names may hold line breaks, and rule_warnings() names them escaped.
            BOOST_LOG_TRIVIAL(info) << "[PreciseSeamRules] nowarn: object_id=" << model_object.id().id
                                    << " volume_id=" << volumes[i]->id().id << " type or order differs from the object list";
    }

    const auto collect = [&](const std::vector<size_t> &group, std::vector<const ModelVolume*> &out) {
        for (size_t i : group)
            if (layers[i]) {
                effective_out.emplace(volumes[i], EffectiveModifier{ruled[i].effective_type(), std::move(*layers[i])});
                out.push_back(volumes[i]);
            }
    };
    collect(order.strong, strong_volumes_out);
    collect(order.weak, weak_volumes_out);
    has_strong_out = !strong_volumes_out.empty();

    // Weak modifiers use last-write-wins, so apply the higher-priority ones last.
    std::reverse(weak_volumes_out.begin(), weak_volumes_out.end());
}

std::vector<std::string> rule_warnings(const std::vector<ModifierDiagnostic> &diagnostics)
{
    // User text is escaped: a line break would hide the rest of the message in the export dialog.
    const auto single_line = [](const std::string &text) {
        std::string out;
        out.reserve(text.size());
        for (const char c : text) {
            if (c == '\n')
                out += "\\n";
            else if (c == '\r')
                out += "\\r";
            else if (c == '\t')
                out += "\\t";
            else if (static_cast<unsigned char>(c) < 0x20 || c == 0x7f)
                out += (boost::format("\\x%02X") % int(static_cast<unsigned char>(c))).str();
            else
                out += c;
        }
        return out;
    };
    // Where an element is, so that equal rules in different sources stay distinct messages.
    const auto location = [](const RuleDiagnostic &d) {
        return d.source.field ? (boost::format(_u8L("rules field, element %1%")) % d.element).str() :
                                (boost::format(_u8L("name block %1%, element %2%")) % d.source.name_block % d.element).str();
    };
    // Messages in order, with the number of diagnostics behind each: modifiers or objects with equal names
    // give equal texts, which are merged with their count instead of being dropped.
    std::vector<std::pair<std::string, size_t>> messages;
    for (size_t i = 0; i < diagnostics.size(); ++i) {
        const ModifierDiagnostic &entry = diagnostics[i];
        const RuleDiagnostic     &d     = entry.diagnostic;
        // Print objects of one model object report the same diagnostic once.
        const bool repeated = std::any_of(diagnostics.begin(), diagnostics.begin() + i, [&entry, &d](const ModifierDiagnostic &other) {
            return other.volume == entry.volume && other.diagnostic.code == d.code && other.diagnostic.source == d.source &&
                   other.diagnostic.element == d.element && other.diagnostic.text == d.text;
        });
        if (repeated)
            continue;
        const std::string volume = single_line(entry.volume->name);
        const std::string object = single_line(entry.volume->get_object()->name);
        const std::string text   = single_line(d.text);
        BOOST_LOG_TRIVIAL(warning) << "[PreciseSeamRules] object=\"" << object << "\" object_id=" << entry.volume->get_object()->id().id
            << " modifier=\"" << volume << "\" volume_id=" << entry.volume->id().id << " code=" << int(d.code)
            << " source=" << (d.source.field ? std::string("field") : "name block " + std::to_string(d.source.name_block))
            << " element=" << d.element << " text=\"" << text << "\"";
        using Code = RuleDiagnosticCode;
        std::string message;
        switch (d.code) {
        case Code::BrokenMarker:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name.
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\" has a broken rule block in its name, the block was ignored.")) % volume % object).str();
            break;
        case Code::MultipleBlocks:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name.
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\" has several rule blocks in its name, they are applied in order.")) % volume % object).str();
            break;
        case Code::InvalidElement:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": rule \"%3%\" (%4%) is not recognized and was ignored.")) % volume % object % text % location(d)).str();
            break;
        case Code::NumberOutOfRange:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2", %5% the largest allowed number.
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": rule \"%3%\" (%4%) has a number above %5% and was ignored.")) % volume % object % text % location(d) % rule_number_limit).str();
            break;
        case Code::TypeConflict:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": type \"%3%\" (%4%) differs from the first type in its rules and was ignored.")) % volume % object % text % location(d)).str();
            break;
        case Code::PriorityConflict:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": priority \"%3%\" (%4%) differs from the first priority in its rules and was ignored.")) % volume % object % text % location(d)).str();
            break;
        case Code::EmptySet:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": rule \"%3%\" (%4%) selects no existing layer.")) % volume % object % text % location(d)).str();
            break;
        case Code::TypeDiffersFromTree:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name.
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\" acts with the type from its rules, not the one in the object list.")) % volume % object).str();
            break;
        case Code::OrderDiffersFromTree:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name.
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\" acts in the order of its rule priority, not its position in the object list.")) % volume % object).str();
            break;
        }
        // Equal texts also come from two blocks broken the same way, or equal rules in a source.
        const auto same = std::find_if(messages.begin(), messages.end(), [&message](const auto &m) { return m.first == message; });
        if (same == messages.end())
            messages.emplace_back(std::move(message), 1);
        else
            ++same->second;
    }
    std::vector<std::string> out;
    out.reserve(messages.size());
    for (const auto &[message, count] : messages)
        // TRN: Precise Seam layer rules. %1% is a rule warning, %2% how many times it occurred.
        out.push_back(count == 1 ? message : (boost::format(_u8L("%1% (%2% in total)")) % message % count).str());
    return out;
}

} // namespace PreciseSeam
} // namespace Slic3r
