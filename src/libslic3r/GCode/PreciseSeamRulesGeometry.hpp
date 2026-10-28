#ifndef slic3r_PreciseSeamRulesGeometry_hpp_
#define slic3r_PreciseSeamRulesGeometry_hpp_

#include <cstddef>
#include <optional>
#include <unordered_map>
#include <vector>
#include "PreciseSeam.hpp"

// Precise Seam with layer rules: the seam point of a perimeter for the modifiers of an object that uses rules.
// Objects without rules keep the drivers of PreciseSeam.hpp. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
namespace PreciseSeam {

// The point of a perimeter that a magnet modifier acts on.
struct MagnetPoint {
    Point             point;
    PerimeterPosition position; // On the unchanged perimeter.
};

// Point of a prepared, unchanged perimeter nearest to the boundaries of the modifier's regions, or nullopt
// when there are none or all are farther than `limit` (scaled units, compared with the snapping tolerance).
// Equally near points go to the one nearest the centroid of its modifier loop, then of the perimeter,
// then to the rear and left one. With `apex`, the point then moves to the vertex of a corner that the modifier
// stands in, when there is one within reach.
std::optional<MagnetPoint> find_magnet_point(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                             std::optional<double> limit, bool apex = false);

// The point of a magnet with apex whose only crossing is shorter than its length, a touch: the point the type
// takes on that crossing (its start for Left, its end for Right, its middle for Center and the weak types),
// then the tip of the corner the modifier stands in as find_magnet_point() looks for it.
std::optional<MagnetPoint> find_touch_point(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                            const PerimeterSegment &touch, ModelVolumeType type, std::optional<double> limit,
                                            bool apex);

// magnet apex across crossings: where every region of the modifier that crosses the perimeter does so two or
// more times, the tip of the corner between the middles of its crossings, on the shortest part of the
// perimeter that holds them all; the tip nearest the centroid of its region wins. nullopt where the type acts:
// fewer than two crossings, a region crossed once, a discarded fragment, or no vertex deeper than the ends.
std::optional<MagnetPoint> find_crossing_apex(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                              const SegmentExtraction &extracted);

// Result of weak modifier segment processing
struct WeakModifierSegment {
    SeamPlacerImpl::EnforcedBlockedSeamPoint type;  // Enforced/Blocked/Neutral
    Point left_point;               // Coordinates of left (first) point of segment
    PerimeterPosition left_position; // Position on the source perimeter before insertion/refinement.
    Point right_point;              // Coordinates of right (last) point of segment
    PerimeterPosition right_position; // Retained provenance, not an index into the modified polygon.
    // Full containment of a weak modifier: the zone is the whole perimeter, without boundaries
    // (the points and positions above are unused and nothing is inserted for it).
    bool whole_perimeter = false;
    const ModelVolume *modifier = nullptr; // Source modifier, for the warnings.
};

// One Precise Seam modifier as it acts on one print object, set by its layer rules (PreciseSeamRulesApply.hpp).
struct EffectiveModifier {
    ModelVolumeType   type;   // Rule label, or the type in the object list.
    std::vector<bool> layers; // Active object layers by index without raft; empty means every layer.
    bool                  magnet = false; // magnet: a perimeter it does not cross gets its nearest point.
    std::optional<double> magnet_limit;   // Scaled units; none with magnet means any distance.
    bool                  magnet_apex = false; // magnet apex: from the nearest point to the vertex of the corner.
    double                magnet_apex_length = 0.; // Scaled: with apex, one crossing shorter than this is a touch.

    bool active_on(size_t layer_id) const { return layers.empty() || (layer_id < layers.size() && layers[layer_id]); }
};
// Modifiers missing here act with their own type on every layer.
using EffectiveModifiers = std::unordered_map<const ModelVolume*, EffectiveModifier>;

// insert_strong_seam_point() for an object with layer rules: each modifier acts with its effective type on its
// active layers, and a magnet acts where it does not cross the perimeter. A `quiet` call is a probe of the relay
// pre-pass: it logs nothing and takes no `warnings`.
std::optional<Point> insert_strong_seam_point_by_rules(const std::vector<const ModelVolume*> &strong_volumes, Polygon &polygon,
                                                       const PreparedPerimeter &prepared, const Layer *layer,
                                                       const ModifierRegionsCache &slices_cache, PreciseSeamWarnings *warnings,
                                                       const EffectiveModifiers &effective, bool quiet = false);

// collect_weak_modifier_segments() for an object with layer rules, as insert_strong_seam_point_by_rules().
std::vector<WeakModifierSegment> collect_weak_modifier_segments_by_rules(const std::vector<const ModelVolume*> &weak_volumes,
                                                                         Polygon &polygon, const PreparedPerimeter &prepared,
                                                                         const Layer *layer, const ModifierRegionsCache &slices_cache,
                                                                         PreciseSeamWarnings *warnings, const EffectiveModifiers &effective);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_PreciseSeamRulesGeometry_hpp_
