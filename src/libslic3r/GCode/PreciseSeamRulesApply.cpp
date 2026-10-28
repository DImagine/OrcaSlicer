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
#include <unordered_map>
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
#include "libslic3r/libslic3r.h"
#include "PreciseSeam.hpp"
#include "PreciseSeamRulesGeometry.hpp"
#include "PreciseSeamRules.hpp"

namespace Slic3r {
namespace PreciseSeam {

bool uses_precise_seam_rules(const ModelObject &object)
{
    return std::any_of(object.volumes.begin(), object.volumes.end(), [](const ModelVolume *volume) {
        return volume->is_precise_seam() && (!volume->precise_seam_rules.empty() || has_rule_marker(volume->name));
    });
}

const std::vector<const ModelVolume*> &RuleContext::strong_on(size_t layer, const std::vector<const ModelVolume*> &strong,
                                                              std::vector<const ModelVolume*> &buffer) const
{
    if (!relay)
        return strong;
    relay_strong_volumes(*relay, effective, layer, relay->entries.size(), buffer);
    return buffer;
}

LayerNumbering object_numbering(const PrintObject &object, const std::unordered_map<const PrintObject*, LayerNumbering> *plan)
{
    if (plan != nullptr) {
        // The export numbers every object that uses rules. Asserted, not checked at run time: an invariant that
        // cannot fail in the export is documented, not handled.
        const auto it = plan->find(&object);
        assert(it != plan->end());
        if (it != plan->end())
            return it->second;
    }
    // Without the export plan (tests) the object numbers its own layers from 1.
    LayerNumbering own;
    own.count = object.layers().size();
    for (size_t i = 0; i < own.count; ++i)
        own.object_layers.push_back(i + 1);
    return own;
}

namespace {

// Relay groups of an object by label, in volume order and without off or blank volumes. A group is all strong or
// all weak; a mixed one acts without relay, and its members keep their other rules.
struct RelayGroups {
    std::vector<std::optional<size_t>> group_of; // Per volume.
    std::vector<bool>                  mixed;    // Per volume: its label names a mixed group.
    std::vector<uint64_t>              seeds;    // Per group; relay-random: from every member, before layers drop any.
};

RelayGroups find_relay_groups(const std::vector<LayerRules> &rules, const std::vector<RuledModifier> &ruled)
{
    RelayGroups groups{std::vector<std::optional<size_t>>(rules.size()), std::vector<bool>(rules.size(), false), {}};
    std::map<RelayLabel, std::vector<size_t>> labelled;
    for (size_t i = 0; i < rules.size(); ++i)
        if (!rules[i].inert() && rules[i].relay)
            labelled[*rules[i].relay].push_back(i);
    for (const auto &[label, members] : labelled) {
        const auto strong = [&ruled](size_t i) { return is_precise_seam_strong(ruled[i].effective_type()); };
        const bool same   = std::all_of(members.begin(), members.end(), strong) || std::none_of(members.begin(), members.end(), strong);
        for (size_t i : members)
            if (same)
                groups.group_of[i] = groups.seeds.size();
            else
                groups.mixed[i] = true;
        if (!same)
            continue;
        std::vector<std::optional<uint32_t>> seeds;
        for (size_t i : members)
            seeds.push_back(rules[i].relay_seed);
        groups.seeds.push_back(label.random ? relay_random_group_seed(label.number, seeds) : 0);
    }
    return groups;
}

// Height bands of the numbered model layers, matched by the original layer index (numbered_layers() drops that
// index); built only when some volume selects by height.
std::vector<LayerBand> height_bands(const std::vector<LayerRules> &rules, const LayerNumbering &numbering,
                                    const ConstLayerPtrsAdaptor &model_layers)
{
    std::vector<LayerBand> bands;
    const bool uses_heights = std::any_of(rules.begin(), rules.end(), [](const LayerRules &r) {
        return std::any_of(r.operations.begin(), r.operations.end(),
                           [](const RuleLayerOperation &op) { return op.set.kind == RuleLayerSet::Kind::Height; });
    });
    if (!uses_heights)
        return bands;
    // Both come from the same print object, so they always have one entry per layer.
    assert(numbering.object_layers.size() == model_layers.size());
    for (size_t layer = 0; layer < numbering.object_layers.size(); ++layer)
        if (numbering.object_layers[layer] != 0)
            bands.push_back({numbering.object_layers[layer], model_layers[layer]->bottom_z(), model_layers[layer]->print_z});
    return bands;
}

// The numbered layers of an object and their heights, as select_layers() takes them.
struct NumberedLayers {
    const LayerNumbering  &numbering;
    std::vector<size_t>    numbered;
    std::vector<LayerBand> bands;
    double                 model_bottom; // Bottom of the first model layer, above the raft and its gap.
};

// Object layers of one volume: an empty row means every layer; nullopt leaves the volume out: switched off or
// blank, or its rules exclude every layer of this object. Diagnostics of the selection go to `diagnostics`.
std::optional<std::vector<bool>> volume_layers(const LayerRules &rules, const NumberedLayers &layers, std::vector<RuleDiagnostic> &diagnostics)
{
    if (rules.inert())
        return std::nullopt;
    if (!rules.restricts_layers())
        return std::vector<bool>();
    const LayerNumbering &numbering = layers.numbering;
    LayerSelection selection = select_layers(rules, numbering.count, layers.numbered, layers.bands, numbering.object_layers, layers.model_bottom);
    diagnostics.insert(diagnostics.end(), std::make_move_iterator(selection.diagnostics.begin()), std::make_move_iterator(selection.diagnostics.end()));
    std::vector<bool> row(numbering.object_layers.size(), false);
    bool              any = false;
    for (size_t layer = 0; layer < row.size(); ++layer) {
        row[layer] = selection.contains(numbering.object_layers[layer]);
        any        = any || row[layer];
    }
    if (!any)
        return std::nullopt;
    return row;
}

bool in_list(const std::vector<size_t> &list, size_t i) { return std::find(list.begin(), list.end(), i) != list.end(); }

// Diagnostics `own` of volume `i` in source order (the name as a whole, its blocks, then the field) and element
// order; parsing and empty-set diagnostics of one element keep their order. Differences from the object list
// concern the volume as a whole and follow them.
void report_volume(const ModelObject &object, size_t i, const ModelVolume *volume, const LayerRules &rules, std::vector<RuleDiagnostic> own,
                   const EffectiveOrder &order, bool mixed, std::vector<ModifierDiagnostic> &out)
{
    const auto source_rank = [](const RuleDiagnostic &d) {
        return std::pair<size_t, size_t>(d.source.field ? std::numeric_limits<size_t>::max() : d.source.name_block, d.element);
    };
    std::stable_sort(own.begin(), own.end(),
                     [&source_rank](const RuleDiagnostic &a, const RuleDiagnostic &b) { return source_rank(a) < source_rank(b); });
    if (in_list(order.type_changed, i))
        own.push_back({RuleDiagnosticCode::TypeDiffersFromTree, {}, 0, {}});
    if (in_list(order.reordered, i))
        own.push_back({RuleDiagnosticCode::OrderDiffersFromTree, {}, 0, {}});
    if (mixed)
        own.push_back({RuleDiagnosticCode::RelayGroupMixed, {}, 0, relay_label_text(*rules.relay)});
    for (RuleDiagnostic &diagnostic : own)
        out.push_back({volume, std::move(diagnostic)});
    if (rules.name_overridden)
        BOOST_LOG_TRIVIAL(info) << "[PreciseSeamRules] override: object_id=" << object.id().id << " volume_id=" << volume->id().id
                                << " the rules field replaces the name rules";
    // nowarn hides the differences from the user, not from the log.
    if (in_list(order.suppressed, i))
        // IDs only: names may hold line breaks, and rule_warnings() names them escaped.
        BOOST_LOG_TRIVIAL(info) << "[PreciseSeamRules] nowarn: object_id=" << object.id().id << " volume_id=" << volume->id().id
                                << " type or order differs from the object list";
}

// The effective modifier of a volume that takes part, on `layers`.
EffectiveModifier effective_modifier(const RuledModifier &ruled, const LayerRules &rules, std::vector<bool> layers)
{
    EffectiveModifier modifier{ruled.effective_type(), std::move(layers)};
    if (rules.magnet) {
        modifier.magnet             = true;
        modifier.magnet_apex        = rules.magnet->apex;
        modifier.magnet_apex_length = scale_(double(rules.magnet->apex_length) / nanometres_per_mm);
        // Kept in scaled units as a double: a limit below one unit must not round to a touch or to none.
        if (rules.magnet->limit)
            modifier.magnet_limit = scale_(double(*rules.magnet->limit) / nanometres_per_mm);
    }
    return modifier;
}

// An effective order as places: a volume on its own, or the block of its group. A group stands at the place of its
// most senior member, even one whose rules exclude every layer. A block keeps the members that take part, in the
// object list order, which is the queue order.
struct Blocks {
    struct Place {
        size_t                volume;
        std::optional<size_t> block;
    };
    std::vector<std::optional<size_t>> block_of; // Per group.
    std::vector<std::vector<size_t>>   members;  // Per block.

    std::vector<Place> place(const std::vector<size_t> &list, const RelayGroups &groups, const std::vector<std::optional<std::vector<bool>>> &layers)
    {
        std::vector<Place> places;
        for (size_t i : list)
            if (!groups.group_of[i]) {
                if (layers[i])
                    places.push_back({i, std::nullopt});
            } else {
                std::optional<size_t> &block = block_of[*groups.group_of[i]];
                if (!block) {
                    block = members.size();
                    members.emplace_back();
                    places.push_back({i, block});
                }
                if (layers[i])
                    members[*block].push_back(i);
            }
        for (const Place &place : places)
            if (place.block)
                std::sort(members[*place.block].begin(), members[*place.block].end());
        return places;
    }
};

// Layers held by each member of a weak group, in full rows, as an empty one means every layer. Weak members take
// the turn without a pre-pass: a layer goes to its holder among the members whose rules keep it, whatever the
// zones meet. Rows without one stay, as relay-random hashes the real index.
std::vector<std::vector<bool>> weak_relay_layers(const std::vector<size_t> &members, const std::vector<std::optional<std::vector<bool>>> &layers,
                                                 size_t layer_count, const RelayLabel &label, uint64_t seed)
{
    std::vector<std::vector<uint8_t>> active(layer_count, std::vector<uint8_t>(members.size()));
    for (size_t layer = 0; layer < active.size(); ++layer)
        for (size_t m = 0; m < members.size(); ++m)
            active[layer][m] = layers[members[m]]->empty() || (*layers[members[m]])[layer];
    const std::vector<size_t> holders = label.random ? relay_random_holders(active, members.size(), seed) :
                                                       relay_holders(active, members.size());
    std::vector<std::vector<bool>> rows(members.size(), std::vector<bool>(holders.size()));
    for (size_t layer = 0; layer < holders.size(); ++layer)
        if (holders[layer] != RelayOrder::no_holder)
            rows[holders[layer]][layer] = true;
    return rows;
}

} // namespace

RuleContext prepare_rule_context(std::vector<const ModelVolume*> &strong_volumes_out,
                                 std::vector<const ModelVolume*> &weak_volumes_out,
                                 bool &has_strong_out,
                                 const ModelObject &model_object,
                                 const LayerNumbering &numbering,
                                 const ConstLayerPtrsAdaptor &model_layers,
                                 std::vector<ModifierDiagnostic> &diagnostics_out)
{
    strong_volumes_out.clear();
    weak_volumes_out.clear();
    has_strong_out = false;
    RuleContext context;

    // Rules are parsed once per volume, before the parallel layer processing.
    std::vector<const ModelVolume*> volumes;
    std::vector<LayerRules>         rules;
    for (const ModelVolume* volume : model_object.volumes)
        if (volume->is_precise_seam()) {
            volumes.push_back(volume);
            rules.push_back(parse_layer_rules(volume->name, volume->precise_seam_rules));
        }
    if (volumes.empty())
        return context;
    std::vector<RuledModifier> ruled;
    for (size_t i = 0; i < volumes.size(); ++i)
        ruled.push_back({volumes[i]->type(), &rules[i]});
    const EffectiveOrder order  = effective_order(ruled);
    const RelayGroups    groups = find_relay_groups(rules, ruled);

    // Layers of each volume, with its diagnostics in volume order. Nothing reads the diagnostics of the rules after.
    const NumberedLayers numbered{numbering, numbering.numbered_layers(), height_bands(rules, numbering, model_layers),
                                  model_layers.empty() ? 0. : model_layers.front()->bottom_z()};
    std::vector<std::optional<std::vector<bool>>> layers(volumes.size());
    for (size_t i = 0; i < volumes.size(); ++i) {
        std::vector<RuleDiagnostic> own = std::move(rules[i].diagnostics);
        layers[i]                       = volume_layers(rules[i], numbered, own);
        report_volume(model_object, i, volumes[i], rules[i], std::move(own), order, groups.mixed[i], diagnostics_out);
    }

    const auto take_part = [&](size_t i, std::vector<const ModelVolume*> &out) {
        context.effective.emplace(volumes[i], effective_modifier(ruled[i], rules[i], std::move(*layers[i])));
        out.push_back(volumes[i]);
    };
    Blocks     blocks{std::vector<std::optional<size_t>>(groups.seeds.size()), {}};
    RelayOrder relay;
    for (const auto &[i, block] : blocks.place(order.strong, groups, layers)) {
        // A single member acts on its own at the place of the group.
        if (!block || blocks.members[*block].size() == 1) {
            const size_t volume = block ? blocks.members[*block].front() : i;
            relay.entries.push_back({volumes[volume]});
            take_part(volume, strong_volumes_out);
        } else if (!blocks.members[*block].empty()) {
            relay.entries.push_back({nullptr, relay.groups.size()});
            RelayOrder::Group &group = relay.groups.emplace_back();
            group.label              = *rules[i].relay; // The anchor is a member, even without layers.
            group.seed               = groups.seeds[*groups.group_of[i]];
            for (size_t member : blocks.members[*block]) {
                group.members.push_back(volumes[member]);
                take_part(member, strong_volumes_out);
            }
        }
    }
    if (!relay.groups.empty())
        context.relay = std::move(relay);
    for (const auto &[i, block] : blocks.place(order.weak, groups, layers)) {
        // A single member acts on its own at the place of the group.
        if (!block || blocks.members[*block].size() == 1) {
            take_part(block ? blocks.members[*block].front() : i, weak_volumes_out);
            continue;
        }
        const std::vector<size_t> &members = blocks.members[*block];
        if (members.empty())
            continue;
        // The anchor is a member, even without layers.
        std::vector<std::vector<bool>> held = weak_relay_layers(members, layers, numbering.object_layers.size(), *rules[i].relay,
                                                                groups.seeds[*groups.group_of[i]]);
        // A member that holds no layer is left out, like one whose rules exclude every layer.
        for (size_t m = 0; m < members.size(); ++m)
            if (std::find(held[m].begin(), held[m].end(), true) != held[m].end()) {
                *layers[members[m]] = std::move(held[m]);
                take_part(members[m], weak_volumes_out);
            }
    }
    has_strong_out = !strong_volumes_out.empty();

    // Weak modifiers use last-write-wins, so apply the higher-priority ones last.
    std::reverse(weak_volumes_out.begin(), weak_volumes_out.end());
    return context;
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

uint64_t relay_random_group_seed(int64_t group_number, const std::vector<std::optional<uint32_t>> &seeds)
{
    uint64_t seed = mix(uint64_t(group_number));
    if (std::none_of(seeds.begin(), seeds.end(), [](const std::optional<uint32_t> &s) { return s.has_value(); }))
        return seed;
    // One step per member keeps its place; +1 in 64 bits tells seed0 from no seed, also for 2^32 - 1.
    for (const std::optional<uint32_t> &s : seeds)
        seed = mix(seed ^ (s ? uint64_t(*s) + 1 : 0));
    return seed;
}

std::vector<size_t> relay_random_holders(const std::vector<std::vector<uint8_t>> &wins, size_t member_count, uint64_t group_seed)
{
    std::vector<size_t> holders(wins.size(), RelayOrder::no_holder);
    size_t              previous = RelayOrder::no_holder; // Holder of the last layer with a winner.
    const uint64_t      seed     = group_seed;
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

void schedule_relay(RuleContext &context, const PrintObject &object, const ModifierRegionsCache &slices,
                    const LayerPerimeters &perimeters, const std::function<void()> &throw_if_canceled)
{
    if (!context.relay)
        return;
    RelayOrder               &order       = *context.relay;
    const EffectiveModifiers &effective   = context.effective;
    const size_t              layer_count = object.layers().size();
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
                        return insert_strong_seam_point_by_rules(volumes, polygon, prepared, &model_layer, slices, nullptr, effective, true)
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
        group.holders = group.label.random ? relay_random_holders(wins, group.members.size(), group.seed) :
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
        case Code::MagnetConflict:
            // TRN: Precise Seam layer rules. %1% is the modifier name, %2% the object name, %3% the rule text, %4% its place, e.g. "rules field, element 2". Do not translate "magnet".
            message = (boost::format(_u8L("Precise Seam: modifier \"%1%\" of \"%2%\": magnet \"%3%\" (%4%) differs from the first magnet in its rules and was ignored.")) % volume % object % text % location(d)).str();
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
