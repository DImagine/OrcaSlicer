#pragma once

#include "PreciseSeam.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include <cstddef>
#include <utility>
#include <vector>

namespace Slic3r::PreciseSeam::detail {

// Binding intermediates retain source edge identity until segment assembly.
struct ClippedEdgeInterval {
    size_t edge;
    double begin;
    double end;
    Point first;
    Point last;
};

struct FragmentBindingFailure {
    size_t pair_index = 0;
    const char *reason = "empty fragment";
};

// Failure rolls back this fragment only; earlier bindings remain intact.
bool append_projected_fragment(const Polyline &fragment, const Polygon &perimeter,
                               std::vector<ClippedEdgeInterval> &intervals,
                               FragmentBindingFailure &failure);

// Exact path, then projection path, without the fallback; intervals are unchanged on failure.
// Exposed so tests can show that a fragment needs the fallback in append_fragment().
bool bind_fragment(const Polyline &fragment, const Polygon &perimeter,
                   std::vector<ClippedEdgeInterval> &intervals, FragmentBindingFailure &failure);

// Binds one fragment. After a failure it tries the rare-case repair and the contact rule, both logged
// as recoveries. Returns false when the fragment is discarded.
bool append_fragment(const Polyline &fragment, const Polygon &perimeter,
                     std::vector<ClippedEdgeInterval> &intervals,
                     const ExtractionContext &context, size_t fragment_index);

} // namespace Slic3r::PreciseSeam::detail

namespace Slic3r::PreciseSeam {

// Shared with the drivers for layer rules (PreciseSeamRulesGeometry.cpp), which insert and refine seam points
// and prepare weak zones the same way.

// Snapping radius for point insertion, shared by every rule that must match insertion. Deliberately in
// scaled units (1 um by default, 10 um on large beds), well above single-precision candidate steps.
static constexpr coord_t TOLERANCE_LINEAR = 1000;

// Inserts a point on edge `edge_start_idx`, or reuses a vertex closer than TOLERANCE_LINEAR.
std::pair<Point, size_t> insert_point_into_perimeter(const Point &point, size_t edge_start_idx, Polygon &perimeter_polygon);
// Adds a helper point TOLERANCE_LINEAR after (+1) or before (-1) vertex `point_idx`.
void refine_at_vertex(size_t point_idx, int direction, Polygon &perimeter_polygon);
// The seam point type of a weak modifier type.
SeamPlacerImpl::EnforcedBlockedSeamPoint convert_weak_modifier_type(ModelVolumeType type);
// Inserts the boundaries of weak zones into the perimeter, in modifier priority order.
std::vector<WeakModifierSegment> prepare_weak_modifier_segments(std::vector<WeakModifierSegment> result, Polygon &polygon);

} // namespace Slic3r::PreciseSeam
