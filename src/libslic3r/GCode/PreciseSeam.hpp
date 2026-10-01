#ifndef slic3r_PreciseSeam_hpp_
#define slic3r_PreciseSeam_hpp_

#include <atomic>
#include <cassert>
#include <optional>
#include <vector>
#include <unordered_map>
#include <utility>
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "SeamPlacer.hpp"

// Both modifier kinds consume ready perimeter segments. Strong chooses the longest
// within the first matching modifier; weak applies all segments in priority order.
// Full containment is skipped by strong and Blocked; Enforced and Neutral then type the whole
// perimeter. Failed fragments are ignored with diagnostics.

namespace Slic3r {
namespace PreciseSeam {

// Import EnforcedBlockedSeamPoint from SeamPlacerImpl namespace for convenience
using SeamPlacerImpl::EnforcedBlockedSeamPoint;

// Geometry and its exterior bounds are prepared together, then treated as read-only.
struct ModifierRegion {
    ExPolygon polygon;
    BoundingBox bounds;

    explicit ModifierRegion(ExPolygon region)
        : polygon(std::move(region)), bounds(polygon.contour.points) {}
};

using ModifierRegions = std::vector<ModifierRegion>;
using ModifierSlices = std::vector<ModifierRegions>;
// Per-volume slices with cached bounds, shared read-only by both modifier kinds.
using ModifierRegionsCache = std::unordered_map<const ModelVolume*, ModifierSlices>;

// Move sliced geometry into the cache without detaching holes or changing layer indices.
ModifierRegions prepare_modifier_regions(ExPolygons regions);
ModifierSlices prepare_modifier_slices(std::vector<ExPolygons> slices);

// Bound diagnostic volume only; every failed or recovered fragment is still counted and handled.
// The same limit applies separately to failure and recovery markers.
inline constexpr size_t failed_fragment_log_limit = 10;

// Shared by all layers and objects in one SeamPlacer::init(); a new pass starts fresh.
struct PreciseSeamWarnings {
    // Masks of the Precise Seam types that caused each warning reason, one bit per type (type_bit()).
    // The user warning lists the types instead of naming modifiers.
    std::atomic<unsigned> multiple_intersections{0}; // Center/Left/Right with several segments on a perimeter.
    std::atomic<unsigned> full_containment{0};       // Skipped for a perimeter fully inside: Center/Left/Right, Blocked.
    std::atomic<unsigned> failed_types{0};           // Types with at least one discarded fragment.
    std::atomic<size_t> failed_fragments{0}; // Total discarded fragments, for the log summary.
    // Fragments saved by the rare-case fallback or accepted as contacts; log only, no user warning.
    // Clipper is deterministic, so a prismatic model can repeat the same case on every layer.
    std::atomic<size_t> recovered_fragments{0};

    // Per-modifier usage for the "had no effect" reason. SeamPlacer::init() registers every modifier
    // before the parallel phase, so workers never change the map itself and only set the flags.
    // Unregistered modifiers (e.g. direct calls in tests) are simply not tracked.
    struct ModifierUsage {
        std::atomic<bool> checked{false}; // Extracted on at least one perimeter.
        std::atomic<bool> reached{false}; // Gave a segment, full containment or a discarded fragment.
    };
    std::unordered_map<const ModelVolume*, ModifierUsage> modifier_usage;

    // Bit of a Precise Seam type in the masks above, in menu order (Center is bit 0).
    static unsigned type_bit(ModelVolumeType type)
    {
        assert(is_precise_seam(type));
        return 1u << (int(type) - int(ModelVolumeType::PRECISE_SEAM_CENTER));
    }
    // Load before fetch_or: most calls find the bit already set, so shared cache lines stay clean.
    static void mark(std::atomic<unsigned> &mask, ModelVolumeType type)
    {
        const unsigned bit = type_bit(type);
        if ((mask.load(std::memory_order_relaxed) & bit) == 0)
            mask.fetch_or(bit, std::memory_order_relaxed);
    }
};

// Optional caller identity for concise diagnostics when an intersection is discarded.
struct ExtractionContext {
    const Layer *layer = nullptr;
    const ModelVolume *modifier = nullptr;
    PreciseSeamWarnings *warnings = nullptr;
};

// Borrows the source polygon; use only until insertion/refinement changes that polygon.
struct PreparedPerimeter {
    const Polygon &polygon;
    BoundingBox bounds;
    Polyline line;
    bool valid = false;

    explicit PreparedPerimeter(const Polygon &perimeter);
};

// A vertex is represented by its outgoing edge and parameter zero, including vertex 0.
struct PerimeterPosition {
    size_t edge_index;
    double parameter;
};

// Prepared against the immutable perimeter, before insertion shifts its edge indices.
struct StrongSeamTarget {
    Point point;
    size_t edge_index;
};

struct PerimeterSegment {
    Polyline polyline;
    // One bound source edge per polyline interval.
    std::vector<size_t> edge_indices;
    PerimeterPosition begin;
    PerimeterPosition end;
    std::optional<StrongSeamTarget> strong_target; // Absent for weak modifiers and full containment.
    // Scaled arc length, calculated only for Center or comparison of multiple strong segments.
    double length = 0.; // Zero means unmeasured for weak, full containment, and a single Left/Right segment.
};

struct SegmentExtraction {
    std::vector<PerimeterSegment> segments;
    bool full_containment = false;
    bool valid = true; // Invalid perimeter input; discarded fragments do not invalidate other segments.
    // Clipped fragments whose binding failed; one segment may consist of several fragments.
    size_t discarded_fragments = 0; // Failed bindings are ignored, with a warning and diagnostic marker.
};

// Clip a prepared, unchanged perimeter against each nearby modifier region.
// Each exterior keeps its holes; disjoint region bounds are rejected before clipping.
// Outer contours and holes use nonzero winding. The perimeter must have at least three
// vertices and no consecutive duplicates; either traversal direction is accepted.
// Non-contained strong segments receive ready mode points; lengths are measured for Center or comparison.
// Weak and fully contained results retain geometry and bindings without preparing strong data.
SegmentExtraction extract_perimeter_segments(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                             ModelVolumeType mode, const ExtractionContext &context = {});

// Result of weak modifier segment processing
struct WeakModifierSegment {
    EnforcedBlockedSeamPoint type;  // Enforced/Blocked/Neutral
    Point left_point;               // Coordinates of left (first) point of segment
    PerimeterPosition left_position; // Position on the source perimeter before insertion/refinement.
    Point right_point;              // Coordinates of right (last) point of segment
    PerimeterPosition right_position; // Retained provenance, not an index into the modified polygon.
    // Full containment of an Enforced or Neutral modifier: the zone is the whole perimeter, without
    // boundaries (the points and positions above are unused and nothing is inserted for it).
    bool whole_perimeter = false;
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
//   prepared          - preparation of this polygon; valid only until it is modified
//   layer             - current layer
//   slices_cache      - pre-sliced modifier polygons (built once in SeamPlacer::init)
// Returns:
//   Coordinates of inserted point (internal units) or std::nullopt if nothing inserted
std::optional<Point> insert_strong_seam_point(
    const std::vector<const ModelVolume*> &strong_volumes,
    Polygon &polygon,
    const PreparedPerimeter &prepared,
    const Layer *layer,
    const ModifierRegionsCache &slices_cache,
    PreciseSeamWarnings* warnings = nullptr);

// Collect all weak modifier segments for a perimeter polygon
// Processes weak modifiers (ENFORCED/BLOCKED/NEUTRAL) and collects segment boundaries
// Passes the ready boundary array to preparation only after all modifiers are collected.
// Also inserts boundary points into the perimeter polygon (sorted by descending source edge and parameter)
// Refines enforced edges by subdividing them into segments ≤ enforcer_oversampling_distance
// Parameters:
//   weak_volumes      - list of weak precise seam modifiers
//   polygon           - perimeter polygon (will be modified with inserted points and refined edges)
//   prepared          - preparation of this polygon, shared with unsuccessful strong processing
//   layer             - current layer
//   slices_cache      - pre-sliced modifier polygons (built once in SeamPlacer::init)
// Returns:
//   Ordered vector of segments with updated coordinates (same order as weak_volumes list)
std::vector<WeakModifierSegment> collect_weak_modifier_segments(
    const std::vector<const ModelVolume*> &weak_volumes,
    Polygon &polygon,
    const PreparedPerimeter &prepared,
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
