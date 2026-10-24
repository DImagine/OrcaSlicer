#ifndef slic3r_PreciseSeamRulesApply_hpp_
#define slic3r_PreciseSeamRulesApply_hpp_

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>
#include "libslic3r/Model.hpp"
#include "libslic3r/Polygon.hpp"
#include "PreciseSeam.hpp"
#include "PreciseSeamRules.hpp"

// Precise Seam layer rules applied to a print: the branch that SeamPlacer and the G-code export take only
// for objects whose Precise Seam volumes use the rules language. Objects without rules keep the code
// path of Precise Seam without rules. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class ConstLayerPtrsAdaptor;
class Layer;
class PrintObject;
namespace PreciseSeam {

// Whether any Precise Seam volume of the object uses the rules language: a non-empty rules field, or a
// name block marker in either direction, in any letter case. Broken or invalid rules count too, so that
// they are reported. The check reads the raw text and parses nothing.
bool uses_precise_seam_rules(const ModelObject &object);

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
    };
    // An ordinary strong volume, or the block of groups[group] when volume is null.
    struct Entry {
        const ModelVolume *volume = nullptr;
        size_t             group  = 0;
    };
    std::vector<Entry> entries; // Effective strong order, highest first.
    std::vector<Group> groups;  // In the order of their blocks.
};

// init_precise_seam_data() for an object that uses layer rules: strong volumes in priority order, weak ones
// in application order, with their effective type and layers in `effective_out` (every Precise Seam volume
// that takes part, also those without rules of their own). Volumes switched off or excluded from every
// layer are left out. Diagnostics are appended to `diagnostics_out` in volume, source and element order.
// `model_layers` are the layers of the print object that `numbering` numbers, in the same order.
// A relay group stands in `strong_volumes_out` at the place of its most senior member, its members in queue
// order. `relay_out` is set, without holders, only when some group has two members or more.
void init_precise_seam_data_by_rules(
    std::vector<const ModelVolume*> &strong_volumes_out,
    std::vector<const ModelVolume*> &weak_volumes_out,
    bool &has_strong_out,
    const ModelObject &model_object,
    const LayerNumbering &numbering,
    const ConstLayerPtrsAdaptor &model_layers,
    EffectiveModifiers &effective_out,
    std::optional<RelayOrder> &relay_out,
    std::vector<ModifierDiagnostic> &diagnostics_out);

// Holders of one group on consecutive object layers from `wins[layer][member]`: non-zero where the member
// places a seam on a perimeter the senior modifiers leave free; an empty row has no winner. The turn starts
// at the first member, the holder is the first winner from the turn on in circular order, and the turn then
// passes to the member after the holder. A layer without a winner keeps the turn.
std::vector<size_t> relay_holders(const std::vector<std::vector<uint8_t>> &wins, size_t member_count);

// Holders of a relay-random group from the same table: a pseudo-random winner other than the holder of the
// last layer with a winner, which holds again only as the only winner. The choice depends on the table, the
// group number and the object layer index alone, so every export and platform gives the same holders.
std::vector<size_t> relay_random_holders(const std::vector<std::vector<uint8_t>> &wins, size_t member_count, int64_t group_number);

// Strong volumes to try on object layer `layer`, from the entries before `end`: a block starts at the holder
// of its group, goes on in queue order and wraps around, without members whose rules exclude the layer;
// a block without a holder is left out. The holders of these groups must be known. Clears `out` first.
void relay_strong_volumes(const RelayOrder &order, const EffectiveModifiers &effective, size_t layer, size_t end,
                          std::vector<const ModelVolume*> &out);

// Outer perimeters of one object layer, as the seam placer gives them to Precise Seam.
using LayerPerimeters = std::function<Polygons(const Layer &layer)>;

// Fills the holders of every group, in block order: a pre-pass, parallel over the layers, finds the members
// that win on each layer against the senior modifiers and the blocks of the senior groups; then the queue
// runs from the bottom layer up. Probes run the strong driver quietly on fresh copies of the perimeters.
void schedule_relay(RelayOrder &order, const PrintObject &object, const EffectiveModifiers &effective,
                    const ModifierRegionsCache &slices, const LayerPerimeters &perimeters,
                    const std::function<void()> &throw_if_canceled);

// Single-line export warnings for the diagnostics of all print objects, in their order. A diagnostic that
// print objects of one model object repeat is reported once; equal texts of different volumes are merged
// with their count.
std::vector<std::string> rule_warnings(const std::vector<ModifierDiagnostic> &diagnostics);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_PreciseSeamRulesApply_hpp_
