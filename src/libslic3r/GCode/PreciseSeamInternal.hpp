#pragma once

#include "PreciseSeam.hpp"

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

// Both binding paths share the same discard and diagnostic policy. After a failure, an end cut
// rounded next to a vertex of the fragment's own chain is snapped to that vertex and binding is
// retried once; a failed fragment shorter than 1 um is accepted as a contact. Both are a safety net
// for rare rounding cases, not part of the normal path, and are logged as recoveries.
bool append_fragment(const Polyline &fragment, const Polygon &perimeter,
                     std::vector<ClippedEdgeInterval> &intervals,
                     const ExtractionContext &context, size_t fragment_index);

} // namespace Slic3r::PreciseSeam::detail
