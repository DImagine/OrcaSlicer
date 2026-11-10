#pragma once

#include <cstddef>
#include <vector>

#include "ExPolygon.hpp"

namespace Slic3r {

class ExtrusionEntityCollection;
class Flow;

// Outer wall printed in several thin passes (option split_outer_wall).
//
// The outer wall of a layer is printed as a stack of thin loops, the passes, before anything else of that layer.
// The inner edge of the stack is vertical and stays where a regular outer wall of the nominal width would end, so the
// inner walls and the infill are generated the regular way and printed at the full layer height afterwards.
// The outer edge of every pass follows the outline of the object at the middle of that pass. The passes therefore
// differ in width and reproduce the slope of the surface in steps of the pass height.
// A layer is printed this way only if all its passes fit the allowed line width range.
struct OuterWallPassSlices
{
    // Outline of the object at the middle of the layer and at the middle of each pass, bottom pass first.
    // The outlines are not XY compensated: only the offset of a pass outline from the middle one is used.
    ExPolygons              mid;
    std::vector<ExPolygons> passes;
    // Allowed width of a pass, unscaled.
    double                  min_width { 0. };
    double                  max_width { 0. };

    size_t count() const { return this->passes.size(); }
    // Width of the outer wall the wall generator works with on such a layer.
    double nominal_width() const { return 0.5 * (this->min_width + this->max_width); }
};

// Whether the passes following the given outlines stay within the allowed width range everywhere.
// With with_holes == false the holes of the outlines are ignored.
bool outer_wall_passes_fit(const ExPolygons              &mid,
                           const std::vector<ExPolygons> &passes,
                           double                         min_width,
                           double                         max_width,
                           bool                           with_holes);

// Replace the external loops of perimeters by their top passes and append the passes below to lower_passes, one
// collection per pass, bottom pass first. pass_flow is the outer wall flow at the pass height.
// Loops which cannot be split (overhangs, another width, holes when with_holes == false) are left as they are.
void split_outer_wall_loops(ExtrusionEntityCollection              &perimeters,
                            std::vector<ExtrusionEntityCollection> &lower_passes,
                            const OuterWallPassSlices              &slices,
                            const Flow                             &pass_flow,
                            bool                                    with_holes);

} // namespace Slic3r
