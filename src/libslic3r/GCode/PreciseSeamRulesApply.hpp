#ifndef slic3r_PreciseSeamRulesApply_hpp_
#define slic3r_PreciseSeamRulesApply_hpp_

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
#include "libslic3r/Model.hpp"
#include "libslic3r/Polygon.hpp"
#include "PreciseSeam.hpp"
#include "PreciseSeamRulesGeometry.hpp"
#include "PreciseSeamRules.hpp"

// Precise Seam layer rules applied to a print: the branch that SeamPlacer and the G-code export take only
// for objects whose Precise Seam volumes use the rules language. Objects without rules keep the code
// path of Precise Seam without rules. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class ConstLayerPtrsAdaptor;
class Layer;
class PrintObject;
namespace PreciseSeam {

// A rule diagnostic tied to its volume, for the export warnings.
struct ModifierDiagnostic {
    const ModelVolume *volume;
    RuleDiagnostic     diagnostic;
};

// Relay groups of one print object: the strong order with each group as one block at the place of its most
// senior member, and the member that holds the turn of each group on each object layer.
struct RelayOrder {
    static constexpr size_t no_holder = size_t(-1);

    struct Group {
        std::vector<const ModelVolume*> members; // Queue order, which is the object list order.
        std::vector<size_t>             holders; // Per object layer: index into members, or no_holder.
        RelayLabel                      label;   // relay-random groups choose their holders differently.
        uint64_t                        seed = 0; // relay-random: relay_random_group_seed() of the group.
    };
    // An ordinary strong volume, or the block of groups[group] when volume is null.
    struct Entry {
        const ModelVolume *volume = nullptr;
        size_t             group  = 0;
    };
    std::vector<Entry> entries; // Effective strong order, highest first.
    std::vector<Group> groups;  // In the order of their blocks.
};

// Layer rules of one print object that uses them: prepared before the seam pass in two stages,
// prepare_rule_context() and schedule_relay(), then only read by the pass, in parallel over the layers.
struct RuleContext {
    // Every Precise Seam volume that takes part, also those without rules of their own.
    EffectiveModifiers effective;
    // Set, with holders after schedule_relay(), only when some strong group has two members or more.
    std::optional<RelayOrder> relay;

    // Strong volumes to try on object layer `layer`: `strong`, the list prepare_rule_context() gave, or with
    // relay groups the order of this layer, built in `buffer`.
    const std::vector<const ModelVolume*> &strong_on(size_t layer, const std::vector<const ModelVolume*> &strong,
                                                     std::vector<const ModelVolume*> &buffer) const;
};

// Rule numbers of `object` in `plan`, the numbering of the export (SeamPlacer::PreciseSeamNumbering), which has
// an entry for every object that uses rules; nullptr without a plan (direct calls, tests).
const LayerNumbering *object_numbering(const PrintObject &object, const std::unordered_map<const PrintObject*, LayerNumbering> *plan);

// init_precise_seam_data() for an object that uses layer rules: strong volumes in priority order, weak ones
// in application order, with their effective type and layers in the context. Volumes switched off or excluded
// from every layer are left out. Diagnostics are appended to `diagnostics_out` in volume, source and element
// order. `model_layers` are the layers of the print object that `plan_numbering` numbers, in the same order;
// without it the object numbers them itself from 1.
// A relay group stands in `strong_volumes_out` at the place of its most senior member, its members in queue
// order; the context gets it without holders. A weak group stands likewise in `weak_volumes_out`, and its
// holders are set here: each member keeps the layers it holds.
RuleContext prepare_rule_context(std::vector<const ModelVolume*> &strong_volumes_out,
                                 std::vector<const ModelVolume*> &weak_volumes_out,
                                 bool &has_strong_out,
                                 const ModelObject &model_object,
                                 const LayerNumbering *plan_numbering,
                                 const ConstLayerPtrsAdaptor &model_layers,
                                 std::vector<ModifierDiagnostic> &diagnostics_out);

// Holders of one group on consecutive object layers from `wins[layer][member]`: non-zero where the member
// places a seam on a perimeter the senior modifiers leave free; an empty row has no winner. The turn starts
// at the first member, the holder is the first winner from the turn on in circular order, and the turn then
// passes to the member after the holder. A layer without a winner keeps the turn.
std::vector<size_t> relay_holders(const std::vector<std::vector<uint8_t>> &wins, size_t member_count);

// Seed of a relay-random group from its number and the seeds of all its members in queue order, those whose
// rules exclude every layer included; `seeds` holds nullopt for a member without one. Without any seed it
// depends on the number alone, as before seeds; otherwise every seed and its place in the queue take part.
uint64_t relay_random_group_seed(int64_t group_number, const std::vector<std::optional<uint32_t>> &seeds);

// Holders of a relay-random group from the same table: a pseudo-random winner other than the holder of the
// last layer with a winner, which holds again only as the only winner. The choice depends on the table, the
// group seed and the object layer index alone, so every export and platform gives the same holders.
std::vector<size_t> relay_random_holders(const std::vector<std::vector<uint8_t>> &wins, size_t member_count, uint64_t group_seed);

// Strong volumes to try on object layer `layer`, from the entries before `end`: a block starts at the holder
// of its group, goes on in queue order and wraps around, without members whose rules exclude the layer;
// a block without a holder is left out. The holders of these groups must be known. Clears `out` first.
void relay_strong_volumes(const RelayOrder &order, const EffectiveModifiers &effective, size_t layer, size_t end,
                          std::vector<const ModelVolume*> &out);

// Outer perimeters of one object layer, as the seam placer gives them to Precise Seam.
using LayerPerimeters = std::function<Polygons(const Layer &layer)>;

// Second stage of the context, once its modifiers are sliced; nothing to do without relay. Fills the holders of
// every group, in block order: a pre-pass, parallel over the layers, finds the members that win on each layer
// against the senior modifiers and the blocks of the senior groups; then the queue runs from the bottom layer
// up. Probes run the strong driver quietly on fresh copies of the perimeters.
void schedule_relay(RuleContext &context, const PrintObject &object, const ModifierRegionsCache &slices,
                    const LayerPerimeters &perimeters, const std::function<void()> &throw_if_canceled);

// Single-line export warnings for the diagnostics of all print objects, in their order. A diagnostic that
// print objects of one model object repeat is reported once; equal texts of different volumes are merged
// with their count.
std::vector<std::string> rule_warnings(const std::vector<ModifierDiagnostic> &diagnostics);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_PreciseSeamRulesApply_hpp_
