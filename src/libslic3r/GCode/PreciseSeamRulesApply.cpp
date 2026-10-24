// Precise Seam layer rules applied to a print. Design: docs/HLSD/precise-seam.md
#include "PreciseSeamRulesApply.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include "libslic3r/I18N.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Polygon.hpp"
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
    std::optional<RelayOrder> &relay_out,
    std::vector<ModifierDiagnostic> &diagnostics_out)
{
    strong_volumes_out.clear();
    weak_volumes_out.clear();
    effective_out.clear();
    relay_out.reset();
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

    // Relay groups by label, in volume order and without off volumes. A group is all strong or all weak; a
    // mixed one acts without relay, and its members keep their other rules.
    std::map<RelayLabel, std::vector<size_t>> labelled;
    for (size_t i = 0; i < volumes.size(); ++i)
        if (!rules[i].off && rules[i].relay)
            labelled[*rules[i].relay].push_back(i);
    std::vector<std::optional<size_t>> group_of(volumes.size());
    std::vector<bool>                  mixed(volumes.size(), false);
    size_t                             group_count = 0;
    for (const auto &[label, members] : labelled) {
        const auto strong = [&ruled](size_t i) { return is_precise_seam_strong(ruled[i].effective_type()); };
        const bool same   = std::all_of(members.begin(), members.end(), strong) || std::none_of(members.begin(), members.end(), strong);
        for (size_t i : members)
            if (same)
                group_of[i] = group_count;
            else
                mixed[i] = true;
        group_count += same ? 1 : 0;
    }

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
        if (mixed[i])
            own.push_back({RuleDiagnosticCode::RelayGroupMixed, {}, 0, relay_label_text(*rules[i].relay)});
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

    const auto take_part = [&](size_t i, std::vector<const ModelVolume*> &out) {
        effective_out.emplace(volumes[i], EffectiveModifier{ruled[i].effective_type(), std::move(*layers[i])});
        out.push_back(volumes[i]);
    };
    // A group stands at the place of its most senior member, even one whose rules exclude every layer.
    // A block keeps the members that take part, in the object list order, which is the queue order.
    std::vector<std::optional<size_t>> block_of(group_count);
    std::vector<std::vector<size_t>>   blocks;
    // An effective order as places: a volume, or the block placed there.
    const auto place = [&](const std::vector<size_t> &list) {
        std::vector<std::pair<size_t, std::optional<size_t>>> places;
        for (size_t i : list)
            if (!group_of[i]) {
                if (layers[i])
                    places.emplace_back(i, std::nullopt);
            } else {
                std::optional<size_t> &block = block_of[*group_of[i]];
                if (!block) {
                    block = blocks.size();
                    blocks.emplace_back();
                    places.emplace_back(i, block);
                }
                if (layers[i])
                    blocks[*block].push_back(i);
            }
        for (const auto &[i, block] : places)
            if (block)
                std::sort(blocks[*block].begin(), blocks[*block].end());
        return places;
    };
    RelayOrder relay;
    for (const auto &[i, block] : place(order.strong)) {
        // A single member acts on its own at the place of the group.
        if (!block || blocks[*block].size() == 1) {
            const size_t volume = block ? blocks[*block].front() : i;
            relay.entries.push_back({volumes[volume]});
            take_part(volume, strong_volumes_out);
        } else if (!blocks[*block].empty()) {
            relay.entries.push_back({nullptr, relay.groups.size()});
            RelayOrder::Group &group = relay.groups.emplace_back();
            group.label              = *rules[i].relay; // The anchor is a member, even without layers.
            for (size_t member : blocks[*block]) {
                group.members.push_back(volumes[member]);
                take_part(member, strong_volumes_out);
            }
        }
    }
    if (!relay.groups.empty())
        relay_out = std::move(relay);
    for (const auto &[i, block] : place(order.weak)) {
        // A single member acts on its own at the place of the group.
        if (!block || blocks[*block].size() == 1) {
            take_part(block ? blocks[*block].front() : i, weak_volumes_out);
            continue;
        }
        const std::vector<size_t> &members = blocks[*block];
        if (members.empty())
            continue;
        // Weak members take the turn without a pre-pass: a layer goes to its holder among the members whose
        // rules keep it, whatever the zones meet. Rows without one stay, as relay-random hashes the real index.
        std::vector<std::vector<uint8_t>> active(numbering.object_layers.size(), std::vector<uint8_t>(members.size()));
        for (size_t layer = 0; layer < active.size(); ++layer)
            for (size_t m = 0; m < members.size(); ++m)
                active[layer][m] = layers[members[m]]->empty() || (*layers[members[m]])[layer];
        const RelayLabel         &label   = *rules[i].relay; // The anchor is a member, even without layers.
        const std::vector<size_t> holders = label.random ? relay_random_holders(active, members.size(), label.number) :
                                                           relay_holders(active, members.size());
        // Each member keeps only the layers it holds, in a full row, as an empty one means every layer. A member
        // that holds none is left out, like one whose rules exclude every layer.
        for (size_t m = 0; m < members.size(); ++m) {
            std::vector<bool> row(holders.size());
            for (size_t layer = 0; layer < row.size(); ++layer)
                row[layer] = holders[layer] == m;
            if (std::find(row.begin(), row.end(), true) != row.end()) {
                *layers[members[m]] = std::move(row);
                take_part(members[m], weak_volumes_out);
            }
        }
    }
    has_strong_out = !strong_volumes_out.empty();

    // Weak modifiers use last-write-wins, so apply the higher-priority ones last.
    std::reverse(weak_volumes_out.begin(), weak_volumes_out.end());
}

std::vector<size_t> relay_holders(const std::vector<std::vector<uint8_t>> &wins, size_t member_count)
{
    std::vector<size_t> holders(wins.size(), RelayOrder::no_holder);
    size_t              turn = 0;
    for (size_t layer = 0; layer < wins.size(); ++layer) {
        const std::vector<uint8_t> &row = wins[layer];
        if (row.empty())
            continue;
        assert(row.size() == member_count);
        for (size_t k = 0; k < member_count; ++k)
            if (const size_t member = (turn + k) % member_count; row[member] != 0) {
                holders[layer] = member;
                turn           = (member + 1) % member_count;
                break;
            }
    }
    return holders;
}

namespace {
// splitmix64 with its increment, on unsigned 64-bit integers only, so that every platform picks the same.
uint64_t mix(uint64_t x)
{
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}
} // namespace

std::vector<size_t> relay_random_holders(const std::vector<std::vector<uint8_t>> &wins, size_t member_count, int64_t group_number)
{
    std::vector<size_t> holders(wins.size(), RelayOrder::no_holder);
    size_t              previous = RelayOrder::no_holder; // Holder of the last layer with a winner.
    const uint64_t      seed     = mix(uint64_t(group_number));
    for (size_t layer = 0; layer < wins.size(); ++layer) {
        const std::vector<uint8_t> &row = wins[layer];
        if (row.empty())
            continue;
        assert(row.size() == member_count);
        size_t winners      = 0;
        bool   previous_won = false;
        for (size_t member = 0; member < member_count; ++member)
            if (row[member] != 0) {
                ++winners;
                previous_won = previous_won || member == previous;
            }
        if (winners == 0)
            continue;
        // The previous holder is left out, unless it is the only winner.
        const bool   skip_previous = previous_won && winners > 1;
        const size_t candidates    = winners - (skip_previous ? 1 : 0);
        size_t       rank          = candidates == 1 ? 0 : size_t(mix(seed ^ uint64_t(layer)) % candidates);
        for (size_t member = 0; member < member_count; ++member)
            if (row[member] != 0 && !(skip_previous && member == previous) && rank-- == 0) {
                holders[layer] = member;
                break;
            }
        previous = holders[layer];
    }
    return holders;
}

void relay_strong_volumes(const RelayOrder &order, const EffectiveModifiers &effective, size_t layer, size_t end,
                          std::vector<const ModelVolume*> &out)
{
    out.clear();
    for (size_t e = 0; e < end; ++e) {
        const RelayOrder::Entry &entry = order.entries[e];
        if (entry.volume != nullptr) {
            out.push_back(entry.volume);
            continue;
        }
        const RelayOrder::Group &group = order.groups[entry.group];
        // Groups are scheduled in block order, and only the blocks before the one being scheduled are read.
        assert(layer < group.holders.size());
        const size_t holder = group.holders[layer];
        if (holder == RelayOrder::no_holder)
            continue;
        for (size_t k = 0; k < group.members.size(); ++k)
            if (const ModelVolume *member = group.members[(holder + k) % group.members.size()]; effective.at(member).active_on(layer))
                out.push_back(member);
    }
}

void schedule_relay(RelayOrder &order, const PrintObject &object, const EffectiveModifiers &effective,
                    const ModifierRegionsCache &slices, const LayerPerimeters &perimeters,
                    const std::function<void()> &throw_if_canceled)
{
    const size_t layer_count = object.layers().size();
    // Groups one after another: the seniors of a group include the blocks of the groups before it.
    for (size_t end = 0; end < order.entries.size(); ++end) {
        if (order.entries[end].volume != nullptr)
            continue;
        RelayOrder::Group                &group = order.groups[order.entries[end].group];
        std::vector<std::vector<uint8_t>> wins(layer_count);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, layer_count), [&](const tbb::blocked_range<size_t> &range) {
            std::vector<const ModelVolume*> seniors;
            std::vector<const ModelVolume*> member(1);
            for (size_t layer = range.begin(); layer < range.end(); ++layer) {
                throw_if_canceled();
                // Members whose rules keep this layer and that have not won on it yet.
                std::vector<size_t> waiting;
                for (size_t m = 0; m < group.members.size(); ++m)
                    if (effective.at(group.members[m]).active_on(layer))
                        waiting.push_back(m);
                if (waiting.empty())
                    continue;
                const Layer &model_layer = *object.layers()[layer];
                relay_strong_volumes(order, effective, layer, end, seniors);
                std::vector<uint8_t> row(group.members.size(), 0);
                for (const Polygon &perimeter : perimeters(model_layer)) {
                    // Each run takes a fresh counterclockwise copy, as the main pass does, since a placed seam
                    // changes the polygon.
                    const auto places = [&](const std::vector<const ModelVolume*> &volumes) {
                        Polygon polygon = perimeter;
                        polygon.make_counter_clockwise();
                        const PreparedPerimeter prepared(polygon);
                        return insert_strong_seam_point(volumes, polygon, prepared, &model_layer, slices, nullptr, &effective, true)
                            .has_value();
                    };
                    // A perimeter that the seniors take is no win for any member.
                    if (!seniors.empty() && places(seniors))
                        continue;
                    for (auto it = waiting.begin(); it != waiting.end();) {
                        member.front() = group.members[*it];
                        if (places(member)) {
                            row[*it] = 1;
                            it       = waiting.erase(it);
                        } else
                            ++it;
                    }
                    if (waiting.empty())
                        break;
                }
                wins[layer] = std::move(row);
            }
        });
        group.holders = group.label.random ? relay_random_holders(wins, group.members.size(), group.label.number) :
                                             relay_holders(wins, group.members.size());
    }
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
        case Code::RelayConflict:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": relay \"%3%\" (%4%) differs from the first relay in its rules and was ignored.")) % volume % object % text % location(d)).str();
            break;
        case Code::RelayGroupMixed:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the relay group, e.g. "relay2". Do not translate "relay".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": group \"%3%\" acts without relay, as it mixes Seam Center, Seam Left or Seam Right modifiers with Seam Enforced, Seam Blocked or Seam Neutral ones.")) % volume % object % text).str();
            break;
        case Code::IdConflict:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2". Do not translate "id".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": id \"%3%\" (%4%) differs from the first id in its rules and was ignored.")) % volume % object % text % location(d)).str();
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
