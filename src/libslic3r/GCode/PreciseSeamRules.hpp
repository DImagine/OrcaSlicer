#ifndef slic3r_PreciseSeamRules_hpp_
#define slic3r_PreciseSeamRules_hpp_

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>
#include "libslic3r/Model.hpp"

// Precise Seam layer rules: a short language in the volume name (`-^ps=...=ps^-` blocks) and in a
// per-volume field that selects the layers, type and priority of one Precise Seam modifier.
// Design: docs/HLSD/precise-seam.md

namespace Slic3r {
namespace PreciseSeam {

// Every numeric literal must not exceed this value; larger ones invalidate their element.
inline constexpr int64_t rule_number_limit = 1000000;

enum class RuleDiagnosticCode {
    BrokenMarker,        // Reopened, stray closing or unclosed name block.
    MultipleBlocks,      // More than one complete name block (informational).
    InvalidElement,      // Element does not match the grammar; ignored.
    NumberOutOfRange,    // Element matches the grammar but a literal exceeds rule_number_limit; ignored.
    TypeConflict,        // A later type label differs from the first one, which stays in effect.
    PriorityConflict,    // A later priority differs from the first one, which stays in effect.
    EmptySet,            // A valid layer set selects no existing layer (informational).
    TypeDiffersFromTree, // Effective type differs from the volume type in the object list.
    OrderDiffersFromTree, // Effective order of the volume differs from its position in the object list.
    RelayConflict,       // A later relay label differs from the first one, which stays in effect.
    RelayGroupNotStrong  // The relay group of the volume has a weak member, so the group acts without relay.
};

// Name blocks are numbered from 1 in name order; the field always comes after them.
struct RuleSource {
    bool   field      = false;
    size_t name_block = 0; // 1-based, unused for the field.

    bool operator==(const RuleSource &rhs) const { return field == rhs.field && name_block == rhs.name_block; }
};

struct RuleDiagnostic {
    RuleDiagnosticCode code;
    RuleSource         source;
    size_t             element = 0; // 1-based element index in its source; 0 when not tied to an element.
    std::string        text;        // Original element text, or the discarded block content for BrokenMarker.
};

// Layer set kept unevaluated: the number of layers is known only when the print is exported.
struct RuleLayerSet {
    enum class Kind { Range, Formula, Odd, Even, Top, All, Height };
    Kind    kind   = Kind::All;
    int64_t first  = 0;     // Range: A; Formula: k (non-zero); Top: K; Height: A in nanometres.
    int64_t last   = 0;     // Range: B; Formula: b; Height: B in nanometres.
    int64_t step   = 1;     // Range only; a single layer N is the range N-N.
    bool    inside = false; // Height only: zin, the whole layer within [A, B]; a single height zA is zA-A.
};

// Heights are kept exactly, in nanometres (1e-6 mm), and converted to mm only when compared.
inline constexpr int64_t nanometres_per_mm = 1000000;

struct RuleLayerOperation {
    bool         exclude = false; // Mode in effect for this set.
    RuleLayerSet set;
    RuleSource   source;
    size_t       element = 0;
    std::string  text;
};

// Parsed rules of one volume, independent of the print.
struct LayerRules {
    std::vector<RuleLayerOperation> operations; // Valid layer sets in application order.
    std::optional<ModelVolumeType>  type;       // First type label over all sources.
    std::optional<int>              priority;   // First priority over all sources: N for prN, -N for pr-N.
    std::optional<int64_t>          relay;      // First relay label over all sources: N for relayN, 0 for relay.
    bool                            off             = false;
    bool                            nowarn          = false;
    bool                            name_overridden = false; // override first in the field: the name took no part.
    std::vector<RuleDiagnostic>     diagnostics; // Parse and conflict diagnostics in source/element order.

    // Rules without any valid layer set act on every layer (compatibility with volumes without rules).
    bool restricts_layers() const { return !operations.empty(); }
};

// Whether the name holds a block marker in either direction, in any letter case, whole block or not.
bool has_rule_marker(std::string_view volume_name);

// Parses the name blocks, then the field, of one Precise Seam volume. When the first element of the field is
// override, the name is not parsed at all.
LayerRules parse_layer_rules(std::string_view volume_name, std::string_view field);

// Height band (bottom_z, print_z] in mm of one numbered model layer of the object.
struct LayerBand {
    size_t number;
    double bottom_z;
    double print_z;
};

struct LayerSelection {
    std::vector<bool>           layers; // Index is the layer number; index 0 is unused and always false.
    std::vector<RuleDiagnostic> diagnostics; // EmptySet only.

    bool contains(size_t layer) const { return layer < layers.size() && layers[layer]; }
};

// Evaluates the rules for layer numbers 1..layer_count. `object_layers` lists, in ascending order, the
// numbers of the layers that hold this object's model layers; topK selects from them. `bands` holds the
// height band of each of these model layers; height sets select from them.
LayerSelection select_layers(const LayerRules &rules, size_t layer_count, const std::vector<size_t> &object_layers,
                             const std::vector<LayerBand> &bands = {});

// Rule numbers of one print object's layers: preview events of the plate when printing by layer, the
// object's own events from 1 when printing by object.
struct LayerNumbering {
    size_t              count = 0;     // Numbers 1..count exist.
    std::vector<size_t> object_layers; // Number of each PrintObject::layers() entry; 0 when it has none.

    // Numbered layers of the object in ascending order, as select_layers() takes them.
    std::vector<size_t> numbered_layers() const;
};

// One Precise Seam volume of an object, in ModelObject::volumes order.
struct RuledModifier {
    ModelVolumeType   actual; // Type in the object list.
    const LayerRules *rules;  // Never null.

    // The rule label wins over the object list.
    ModelVolumeType effective_type() const { return rules->type.value_or(actual); }
};

// Indices into the input; `off` volumes take no part.
struct EffectiveOrder {
    std::vector<size_t> strong;       // Priority order, highest first.
    std::vector<size_t> weak;         // Priority order, highest first; applied in reverse.
    std::vector<size_t> type_changed; // Effective type differs from the object list (without nowarn).
    std::vector<size_t> reordered;    // Their prN changed the order within the group (without nowarn).
    std::vector<size_t> suppressed;   // nowarn volumes that would be in one of the two lists above.
};

// Groups the volumes by effective type, then orders each group pr1, pr2, ..., without priority,
// pr-1, pr-2, ..., ties by the input index. A strong group always precedes a weak one.
EffectiveOrder effective_order(const std::vector<RuledModifier> &modifiers);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_PreciseSeamRules_hpp_
