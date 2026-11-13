#pragma once

#include "PreciseSeam.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/libslic3r.h"
#include <cstddef>
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

// Source of one edge of the surface line: a source edge shifted towards air, or the arc (convex) or
// bridge (concave) around a source vertex. An arc maps to its vertex; a bridge is never real surface.
struct SurfaceElement {
    size_t source; // Source edge or vertex index.
    bool vertex;
};

// Where coverage of the surface ends on one side of a dead stretch of the surface line (around a concave
// corner) and resumes on the other. The two points meet in space, so coverage of both covers the stretch.
struct SurfaceGap {
    size_t end_element;
    Point end;
    size_t start_element;
    Point start;
};

// The reconstructed surface of one perimeter: the closed surface line with the source of every edge,
// and the parts of it that are real surface, at least the shift away from every part of the centerline.
struct SurfaceLine {
    Polygon polygon;                      // No consecutive or closing duplicates.
    std::vector<SurfaceElement> elements; // One per polygon edge.
    std::vector<Point> shifts;            // Integer shift of every source edge, so shifted edges keep parameters.
    Polylines pieces;                     // Real surface, clipped from the closed line.
    std::vector<SurfaceGap> gaps;
    BoundingBox bounds;                   // Of the pieces.
};

// Builds the surface line of a perimeter accepted by PreparedPerimeter (either direction), shift > 0.
SurfaceLine build_surface_line(const Polygon &perimeter, coord_t shift, bool hole);

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
