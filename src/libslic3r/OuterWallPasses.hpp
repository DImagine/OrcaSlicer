#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "ExPolygon.hpp"
#include "Flow.hpp"

namespace Slic3r {

class ExtrusionEntityCollection;

// Outer wall printed in several thin passes (option split_outer_wall). Design: docs/HLSD/outer-wall-passes.md
//
// The outer wall of a layer is printed as a stack of thin loops, the passes, before anything else of that layer.
// The inner edge of the stack is vertical and stays where a regular outer wall of the nominal width would end, so the
// inner walls and the infill are generated the regular way and printed at the full layer height afterwards.
// The outer edge of every pass follows the outline of the object at the middle of that pass. The passes therefore
// differ in width and reproduce the slope of the surface in steps of the pass height.
// A layer is printed this way only if all its passes fit the allowed line width range.
// Where one wall does not suffice, the inner walls behind it are printed in passes as well and share the offset of the
// outline with it: the vertical edge is then the inner edge of the last of these walls.
struct OuterWallPassSlices
{
    // Outline of the object at the middle of the layer and at the middle of each pass, bottom pass first.
    // The outlines are not XY compensated: only the offset of a pass outline from the middle one is used.
    ExPolygons              mid;
    std::vector<ExPolygons> passes;
    // Width of a pass the walls keep to as long as they can, unscaled.
    double                  min_width { 0. };
    double                  max_width { 0. };
    // Width of a pass the walls go up to once all of them reached the range above. Zero stands for that range.
    double                  limit_min_width { 0. };
    double                  limit_max_width { 0. };
    // Number of the walls printed in passes on this layer, counted from the outer one.
    size_t                  walls { 1 };

    size_t count() const { return this->passes.size(); }
    // Width of the outer wall the wall generator works with on such a layer.
    double nominal_width() const { return 0.5 * (this->min_width + this->max_width); }
    double lowest_width() const { return this->limit_min_width > 0. ? std::min(this->limit_min_width, this->min_width) : this->min_width; }
    double highest_width() const { return std::max(this->limit_max_width, this->max_width); }
};

// Inner walls of a layer region, which may be printed in passes behind its outer wall.
struct OuterWallPassInnerWalls
{
    // Flow of an inner wall at the height of a pass.
    Flow   pass_flow;
    // Distance between the middle lines of the outer wall and of the first inner wall, and of two inner walls, unscaled.
    double first_spacing { 0. };
    double spacing { 0. };
};

// How far the outline of a pass may lie inside and outside of the outline at the middle of the layer for the walls
// printed in passes to follow it, unscaled.
struct OuterWallPassReach
{
    double inwards { 0. };
    double outwards { 0. };
};

// Whether the passes following the given outlines stay within the allowed width range everywhere.
// With with_holes == false the holes of the outlines are ignored.
bool outer_wall_passes_fit(const ExPolygons              &mid,
                           const std::vector<ExPolygons> &passes,
                           double                         min_width,
                           double                         max_width,
                           bool                           with_holes);
// The same for the outlines of the passes within the given reach of the middle outline.
bool outer_wall_passes_fit(const ExPolygons              &mid,
                           const std::vector<ExPolygons> &passes,
                           const OuterWallPassReach      &reach,
                           bool                           with_holes);
// Whether no pass steps back from the pass below it by more than max_step, unscaled, which would uncover the joint of
// two walls of that pass. above: what lies over the top pass, if known; it is checked against that pass as well.
bool outer_wall_pass_steps_fit(const std::vector<ExPolygons> &passes, const ExPolygons *above, double max_step, bool with_holes);

// Replace the external loops of perimeters by their top passes and append the passes below to lower_passes, one
// collection per pass, bottom pass first. pass_flow is the outer wall flow at the pass height.
// layer_outline is the outline of the layer the walls were generated in, with the XY compensations: the outer wall
// follows the surface of the object where it runs along it, and the border of another region of the layer elsewhere.
// Loops which cannot be split (another width, holes when with_holes == false) are left as they are. The passes keep
// the role of an overhanging wall where the loop has it, for its speed and cooling, with the flow of a pass.
// With slices.walls above one the loops of the inner walls behind the outer wall are split as well, which takes
// inner_walls.
void split_outer_wall_loops(ExtrusionEntityCollection              &perimeters,
                            std::vector<ExtrusionEntityCollection> &lower_passes,
                            const OuterWallPassSlices              &slices,
                            const ExPolygons                       &layer_outline,
                            const Flow                             &pass_flow,
                            bool                                    with_holes,
                            const OuterWallPassInnerWalls          *inner_walls = nullptr,
                            // TEMP diagnostics: the loops of the outer wall left as they are, and why.
                            std::vector<std::string>               *rejections = nullptr);

} // namespace Slic3r
