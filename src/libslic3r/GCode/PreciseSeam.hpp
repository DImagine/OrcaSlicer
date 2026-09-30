#ifndef slic3r_PreciseSeam_hpp_
#define slic3r_PreciseSeam_hpp_

#include <atomic>
#include <optional>
#include <vector>
#include <unordered_map>
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "SeamPlacer.hpp"

// Both modifier kinds consume ready perimeter segments. Strong chooses the longest
// within the first matching modifier; weak applies all segments in priority order.
// Full containment is skipped; failed fragments are ignored with diagnostics.

namespace Slic3r {
namespace PreciseSeam {

// Import EnforcedBlockedSeamPoint from SeamPlacerImpl namespace for convenience
using SeamPlacerImpl::EnforcedBlockedSeamPoint;

// Per-volume structured slices, built once and shared read-only by both modifier kinds.
using ModifierRegionsCache = std::unordered_map<const ModelVolume*, std::vector<ExPolygons>>;

// Warning flags set during Precise Seam processing (thread-safe)
struct PreciseSeamWarnings {
    std::atomic<bool> multiple_intersections{false};  // modifier intersects perimeter in multiple separate places (strong only)
    std::atomic<bool> full_containment{false};        // modifier fully contains perimeter, no intersection edges
    std::atomic<bool> intersection_processing_failed{false}; // at least one unbindable fragment was ignored
};

// Optional caller identity for concise diagnostics when an intersection is discarded.
struct ExtractionContext {
    const Layer *layer = nullptr;
    const ModelVolume *modifier = nullptr;
    PreciseSeamWarnings *warnings = nullptr;
};

// A vertex is represented by its outgoing edge and parameter zero, including vertex 0.
struct PerimeterPosition {
    size_t edge_index;
    double parameter;
};

struct PerimeterSegment {
    Polyline polyline;
    // One bound source edge per polyline interval.
    std::vector<size_t> edge_indices;
    PerimeterPosition begin;
    PerimeterPosition end;
    std::vector<double> edge_lengths; // Reused by Center; one length per polyline interval.
    double length = 0.; // Euclidean arc length in scaled coordinates, not squared length.
};

struct SegmentExtraction {
    std::vector<PerimeterSegment> segments;
    bool full_containment = false;
    bool valid = true; // Invalid perimeter input; discarded fragments do not invalidate other segments.
    size_t discarded_segments = 0; // Failed bindings are ignored, with a warning and diagnostic marker.
};

// Clip an immutable, implicitly closed perimeter against each nearby modifier region.
// Each exterior keeps its holes; disjoint region bounds are rejected before clipping.
// Outer contours and holes use nonzero winding. The perimeter must have at least three
// vertices and no consecutive duplicates; either traversal direction is accepted.
SegmentExtraction extract_perimeter_segments(const Polygon &perimeter, const ExPolygons &modifier,
                                             const ExtractionContext &context = {});

// Result of weak modifier segment processing
struct WeakModifierSegment {
    EnforcedBlockedSeamPoint type;  // Enforced/Blocked/Neutral
    Point left_point;               // Coordinates of left (first) point of segment
    size_t left_idx;                // Perimeter vertex index for left_point
    Point right_point;              // Coordinates of right (last) point of segment
    size_t right_idx;               // Perimeter vertex index for right_point
};

// Initialize Precise Seam data by populating provided vectors and flag
// Collects precise seam modifiers and fills output parameters
// Call once during SeamPlacer::init() before gather_seam_candidates()
// Parameters:
//   strong_volumes_out - output vector for strong modifiers (CENTER/LEFT/RIGHT)
//   weak_volumes_out   - output vector for weak modifiers (ENFORCED/BLOCKED/NEUTRAL)
//   has_strong_out     - output flag indicating presence of strong modifiers
//   model_object       - model object containing volumes
void init_precise_seam_data(
    std::vector<const ModelVolume*>& strong_volumes_out,
    std::vector<const ModelVolume*>& weak_volumes_out,
    bool& has_strong_out,
    const ModelObject* model_object);

// Insert strong seam point into perimeter polygon
// Processes strong modifiers (CENTER/LEFT/RIGHT) and inserts seam point into polygon
// Parameters:
//   strong_volumes    - list of strong precise seam modifiers
//   polygon           - perimeter polygon (will be modified if point inserted)
//   layer             - current layer
//   slices_cache      - pre-sliced modifier polygons (built once in SeamPlacer::init)
// Returns:
//   Coordinates of inserted point (internal units) or std::nullopt if nothing inserted
std::optional<Point> insert_strong_seam_point(
    const std::vector<const ModelVolume*> &strong_volumes,
    Polygon &polygon,
    const Layer *layer,
    const ModifierRegionsCache &slices_cache,
    PreciseSeamWarnings* warnings = nullptr);

// Collect all weak modifier segments for a perimeter polygon
// Processes weak modifiers (ENFORCED/BLOCKED/NEUTRAL) and collects segment boundaries
// Passes the ready boundary array to preparation only after all modifiers are collected.
// Also inserts boundary points into the perimeter polygon (sorted by descending arc length)
// Refines enforced edges by subdividing them into segments ≤ enforcer_oversampling_distance
// Parameters:
//   weak_volumes      - list of weak precise seam modifiers
//   polygon           - perimeter polygon (will be modified with inserted points and refined edges)
//   layer             - current layer
//   slices_cache      - pre-sliced modifier polygons (built once in SeamPlacer::init)
// Returns:
//   Ordered vector of segments with updated coordinates (same order as weak_volumes list)
std::vector<WeakModifierSegment> collect_weak_modifier_segments(
    const std::vector<const ModelVolume*> &weak_volumes,
    Polygon &polygon,
    const Layer *layer,
    const ModifierRegionsCache &slices_cache,
    PreciseSeamWarnings* warnings = nullptr);

// Apply weak modifier types to perimeter points based on segment boundaries
// Finds boundary points in refined polygon and sets types for points within segments
// Parameters:
//   weak_segments       - segments with boundary coordinates and types
//   result              - layer seams data to modify
//   perimeter           - perimeter info (start/end indices)
//   some_point_enforced - flag to update if Enforced points are set
void apply_weak_modifiers_to_perimeter(
    const std::vector<WeakModifierSegment> &weak_segments,
    PrintObjectSeamData::LayerSeams &result,
    const SeamPlacerImpl::Perimeter &perimeter,
    bool &some_point_enforced);

// Restore precise seam positions that may have been modified by alignment
// Iterates through all perimeters and restores precise_seam_point positions
void restore_precise_seam_positions(std::vector<PrintObjectSeamData::LayerSeams> &layers);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_PreciseSeam_hpp_
