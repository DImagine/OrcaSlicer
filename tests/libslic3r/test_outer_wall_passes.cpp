#include <catch2/catch_all.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/OuterWallPasses.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/libslic3r.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

namespace {

// Square centered at the origin, optionally with a centered square hole.
ExPolygon square(double side, double hole_side = 0.)
{
    auto ring = [](double s) {
        const coord_t half = scaled<coord_t>(0.5 * s);
        return Polygon({ Point(-half, -half), Point(half, -half), Point(half, half), Point(-half, half) });
    };
    ExPolygon out(ring(side));
    if (hole_side > 0.) {
        Polygon hole = ring(hole_side);
        hole.reverse();
        out.holes.emplace_back(std::move(hole));
    }
    return out;
}

// Passes may be 0.4 to 0.6 mm wide, so an outline may differ from the middle one by 0.1 mm.
const double min_width = 0.4;
const double max_width = 0.6;

// Outer wall of the nominal width as one closed loop through the given points.
ExtrusionEntityCollection wall_loop(Points points, bool hole = false)
{
    points.emplace_back(points.front());
    ExtrusionPath path(erExternalPerimeter, 0.1, float(0.5 * (min_width + max_width)), 0.2f);
    path.polyline = Polyline3(Polyline(std::move(points)));
    ExtrusionLoop loop(path, hole ? elrHole : elrDefault);
    loop.inset_idx = 0;
    ExtrusionEntityCollection out;
    out.append(loop);
    return out;
}

// The same along a square of the given side: counter-clockwise around a contour, clockwise around a hole.
ExtrusionEntityCollection wall_loop(double side, bool hole = false)
{
    const coord_t half = scaled<coord_t>(0.5 * side);
    Points        points { Point(-half, -half), Point(half, -half), Point(half, half), Point(-half, half) };
    if (hole)
        std::reverse(points.begin(), points.end());
    return wall_loop(std::move(points), hole);
}

// Inner wall 0.45 mm wide as one loop along a square of the given side: counter-clockwise around a contour,
// clockwise around a hole.
ExtrusionLoop inner_wall_loop(double side, bool hole = false, int inset = 1)
{
    const coord_t half = scaled<coord_t>(0.5 * side);
    Points        points { Point(-half, -half), Point(half, -half), Point(half, half), Point(-half, half) };
    if (hole)
        std::reverse(points.begin(), points.end());
    points.emplace_back(points.front());
    ExtrusionPath path(erPerimeter, 0.1, 0.45f, 0.2f);
    path.polyline = Polyline3(Polyline(std::move(points)));
    ExtrusionLoop loop(path, hole ? elrHole : elrDefault);
    loop.inset_idx = inset;
    return loop;
}

// Half the extent of a loop along X and the width of its widest path.
std::pair<double, double> reach_and_width(const ExtrusionEntity *entity)
{
    const ExtrusionLoop &loop  = *static_cast<const ExtrusionLoop*>(entity);
    const Polygon        polygon = loop.polygon();
    double               reach = 0.;
    double               width = 0.;
    for (const Point &point : polygon.points)
        reach = std::max(reach, std::abs(unscaled<double>(point.x())));
    for (const ExtrusionPath &path : loop.paths)
        width = std::max(width, double(path.width));
    return { reach, width };
}

} // namespace

TEST_CASE("Passes within the width range fit", "[OuterWallPasses]")
{
    const ExPolygons mid { square(20.) };
    // The sides of the squares move by half of the difference of the side lengths.
    REQUIRE(outer_wall_passes_fit(mid, { { square(20.) }, { square(20.) } }, min_width, max_width, true));
    REQUIRE(outer_wall_passes_fit(mid, { { square(20.16) }, { square(19.84) } }, min_width, max_width, true));
}

TEST_CASE("A pass wider or narrower than allowed does not fit", "[OuterWallPasses]")
{
    const ExPolygons mid { square(20.) };
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.3) }, { square(20.) } }, min_width, max_width, true));
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.) }, { square(19.7) } }, min_width, max_width, true));
}

TEST_CASE("The width range is checked with a tolerance of a few microns", "[OuterWallPasses]")
{
    const ExPolygons mid { square(20.) };
    // Exactly at the limit and 2 microns beyond it.
    REQUIRE(outer_wall_passes_fit(mid, { { square(20.2) }, { square(19.8) } }, min_width, max_width, true));
    REQUIRE(outer_wall_passes_fit(mid, { { square(20.204) }, { square(19.796) } }, min_width, max_width, true));
    // 7 microns beyond the limit.
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.214) }, { square(20.) } }, min_width, max_width, true));
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.) }, { square(19.786) } }, min_width, max_width, true));
}

TEST_CASE("A part appearing or vanishing between the passes does not fit", "[OuterWallPasses]")
{
    const ExPolygons mid { square(20.) };
    ExPolygon island = square(2.);
    island.translate(Point(scaled<coord_t>(30.), 0));
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.), island }, { square(20.) } }, min_width, max_width, true));
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.) }, {} }, min_width, max_width, true));
}

TEST_CASE("A loop is split only along the edge of the outlines it follows", "[OuterWallPasses]")
{
    // Number of the low passes made of a wall loop running along a square of 19.5 mm, for the outline of the layer
    // the wall was generated in and the outlines of the model at the middle of the layer and of its two passes.
    auto low_passes = [](const ExPolygon &layer, const ExPolygon &mid, const ExPolygon &low, const ExPolygon &high) {
        OuterWallPassSlices slices;
        slices.mid       = { mid };
        slices.passes    = { { low }, { high } };
        slices.min_width = min_width;
        slices.max_width = max_width;
        ExtrusionEntityCollection              perimeters = wall_loop(19.5);
        std::vector<ExtrusionEntityCollection> lower_passes;
        split_outer_wall_loops(perimeters, lower_passes, slices, { layer }, Flow(0.5f, 0.1f, 0.4f), true);
        return lower_passes.empty() ? size_t(0) : lower_passes.front().entities.size();
    };
    // A solid square of 20 mm narrowing upwards: the loop runs half its width inside the outline.
    REQUIRE(low_passes(square(20.), square(20.), square(20.05), square(19.95)) == 1);
    // A frame 0.2 mm thick in the model, its wall thickened inwards to 1.2 mm by an XY compensation of its hole.
    // The loop runs in the hole of the outlines of the model, and their nearest edge is not the one it follows.
    REQUIRE(low_passes(square(20., 17.6), square(20., 19.6), square(20.05, 19.65), square(19.95, 19.55)) == 0);
    // A frame 0.55 mm thick in the model, thickened the same way, both its edges moving by 0.1 mm per pass. The
    // loop runs inside the frame, nearer to its outer edge at the middle of the layer and at the high pass, but
    // nearer to its hole at the low pass. The widths the outlines ask for are all allowed.
    REQUIRE(outer_wall_passes_fit({ square(20., 18.9) }, { { square(20.2, 19.1) }, { square(19.8, 18.7) } }, min_width, max_width, true));
    REQUIRE(low_passes(square(20., 16.9), square(20., 18.9), square(20.2, 19.1), square(19.8, 18.7)) == 0);
}

TEST_CASE("A loop along the border of another region is split, but not offset there", "[OuterWallPasses]")
{
    // A slab from -10 to 30 mm along X, 20 mm wide, narrowing upwards by 0.025 mm a side per pass. The wall loop of
    // its left part runs along a square of 19.5 mm: three of its sides follow the surface, the right one at
    // 9.75 mm the border of the other part.
    auto slab = [](double grow) {
        const coord_t left  = scaled<coord_t>(-10. - grow);
        const coord_t right = scaled<coord_t>(30. + grow);
        const coord_t half  = scaled<coord_t>(10. + grow);
        return ExPolygon(Polygon({ Point(left, -half), Point(right, -half), Point(right, half), Point(left, half) }));
    };
    OuterWallPassSlices slices;
    slices.mid       = { slab(0.) };
    slices.passes    = { { slab(0.025) }, { slab(-0.025) } };
    slices.min_width = min_width;
    slices.max_width = max_width;
    ExtrusionEntityCollection              perimeters = wall_loop(19.5);
    std::vector<ExtrusionEntityCollection> lower_passes;
    split_outer_wall_loops(perimeters, lower_passes, slices, { slab(0.) }, Flow(0.5f, 0.1f, 0.4f), true);
    REQUIRE(lower_passes.size() == 1);
    REQUIRE(lower_passes.front().entities.size() == 1);
    // The middle of the low pass, 0.025 mm wider, lies 0.0125 mm further out along the surface and stays where it
    // was along the border, away from the corners.
    const Polygon polygon = static_cast<const ExtrusionLoop*>(lower_passes.front().entities.front())->polygon();
    double        surface_x = 0.;
    double        border_x  = 0.;
    for (const Point &point : polygon.points) {
        surface_x = std::min(surface_x, unscaled<double>(point.x()));
        if (std::abs(unscaled<double>(point.y())) < 5.)
            border_x = std::max(border_x, unscaled<double>(point.x()));
    }
    REQUIRE_THAT(surface_x, WithinAbs(-9.7625, 0.003));
    REQUIRE_THAT(border_x, WithinAbs(9.75, 0.003));
}

TEST_CASE("The corners of the wall of a hole follow the hole", "[OuterWallPasses]")
{
    // A square hole of 8 mm widening upwards by 0.25 mm a side per pass, in a square of 20 mm. The wall of the
    // hole runs along a square of 8.5 mm, its corners further from the corners of the hole than its sides from
    // the sides of the hole.
    OuterWallPassSlices slices;
    slices.mid       = { square(20., 8.) };
    slices.passes    = { { square(20., 7.5) }, { square(20., 8.5) } };
    slices.min_width = 0.2;
    slices.max_width = 0.8;
    REQUIRE(outer_wall_passes_fit(slices.mid, slices.passes, slices.min_width, slices.max_width, true));
    ExtrusionEntityCollection              perimeters = wall_loop(8.5, true);
    std::vector<ExtrusionEntityCollection> lower_passes;
    split_outer_wall_loops(perimeters, lower_passes, slices, { square(20., 8.) }, Flow(0.5f, 0.1f, 0.4f), true);
    REQUIRE(lower_passes.size() == 1);
    REQUIRE(lower_passes.front().entities.size() == 1);
    // Along the sides the low pass is 0.75 mm wide, its middle 0.125 mm nearer to the middle of the hole. Its
    // corners move with the sides: left where they were, they would reach 4.25 mm.
    const Polygon polygon = static_cast<const ExtrusionLoop*>(lower_passes.front().entities.front())->polygon();
    double        reach   = 0.;
    for (const Point &point : polygon.points)
        reach = std::max(reach, std::max(std::abs(unscaled<double>(point.x())), std::abs(unscaled<double>(point.y()))));
    REQUIRE_THAT(reach, WithinAbs(4.125, 0.003));
}

TEST_CASE("A concave corner of the wall of a contour follows the contour", "[OuterWallPasses]")
{
    // An L of two bars 10 mm wide and 20 mm long, narrowing upwards by 0.05 mm a side per pass. Its outline has
    // one concave corner at (10, 10), where the wall loop is further from the outline than along its sides.
    auto l_shape = [](double inset) {
        const coord_t low  = scaled<coord_t>(inset);
        const coord_t mid  = scaled<coord_t>(10. - inset);
        const coord_t high = scaled<coord_t>(20. - inset);
        return Points{ Point(low, low), Point(high, low), Point(high, mid), Point(mid, mid), Point(mid, high), Point(low, high) };
    };
    OuterWallPassSlices slices;
    slices.mid       = { ExPolygon(Polygon(l_shape(0.))) };
    slices.passes    = { { ExPolygon(Polygon(l_shape(-0.05))) }, { ExPolygon(Polygon(l_shape(0.05))) } };
    slices.min_width = min_width;
    slices.max_width = max_width;
    ExtrusionEntityCollection              perimeters = wall_loop(l_shape(0.25));
    std::vector<ExtrusionEntityCollection> lower_passes;
    split_outer_wall_loops(perimeters, lower_passes, slices, slices.mid, Flow(0.5f, 0.1f, 0.4f), true);
    REQUIRE(lower_passes.size() == 1);
    REQUIRE(lower_passes.front().entities.size() == 1);
    // The middle of the low pass, 0.05 mm wider, lies 0.025 mm further out along the sides of the corner, at
    // 9.775 mm. Its corner moves with them, a little further: left where it was, it would stay at 9.75 mm.
    const Polygon polygon = static_cast<const ExtrusionLoop*>(lower_passes.front().entities.front())->polygon();
    double        corner  = 20.;
    for (const Point &point : polygon.points) {
        const double x = unscaled<double>(point.x());
        const double y = unscaled<double>(point.y());
        if (x > 5. && x < 15. && y > 5. && y < 15.)
            corner = std::min(corner, std::min(x, y));
    }
    REQUIRE_THAT(corner, WithinAbs(9.78, 0.01));
}

TEST_CASE("An inner wall behind the outer wall shares its offset", "[OuterWallPasses]")
{
    // A square of 20 mm narrowing upwards by 0.13 mm a side per pass: more than the 0.1 mm the outer wall takes
    // within 0.4 to 0.6 mm. An inner wall 0.45 mm wide, its middle 0.45 mm behind the middle of the outer wall,
    // takes 0.15 mm outwards and 0.05 mm inwards.
    OuterWallPassSlices slices;
    slices.mid       = { square(20.) };
    slices.passes    = { { square(20.26) }, { square(19.74) } };
    slices.min_width = min_width;
    slices.max_width = max_width;
    slices.walls     = 2;
    const OuterWallPassInnerWalls inner_walls { Flow(0.45f, 0.1f, 0.4f), 0.45, 0.4 };

    SECTION("both walls are printed in passes, each within the width range") {
        ExtrusionEntityCollection perimeters = wall_loop(19.5);
        perimeters.append(inner_wall_loop(18.6));
        std::vector<ExtrusionEntityCollection> lower_passes;
        split_outer_wall_loops(perimeters, lower_passes, slices, slices.mid, Flow(0.5f, 0.1f, 0.4f), true, &inner_walls);
        REQUIRE(lower_passes.size() == 1);
        REQUIRE(lower_passes.front().entities.size() == 2);
        REQUIRE(perimeters.entities.size() == 2);
        // Outwards the walls share 0.13 mm as 0.1 to 0.15, inwards as 0.1 to 0.05. The outer wall moves by what the
        // inner wall takes and by a half of its own share, the inner wall by a half of its share.
        const auto [low_outer_reach, low_outer_width] = reach_and_width(lower_passes.front().entities[0]);
        const auto [low_inner_reach, low_inner_width] = reach_and_width(lower_passes.front().entities[1]);
        const auto [top_outer_reach, top_outer_width] = reach_and_width(perimeters.entities[0]);
        const auto [top_inner_reach, top_inner_width] = reach_and_width(perimeters.entities[1]);
        REQUIRE_THAT(low_outer_width, WithinAbs(0.5 + 0.052, 0.003));
        REQUIRE_THAT(low_inner_width, WithinAbs(0.45 + 0.078, 0.003));
        REQUIRE_THAT(top_outer_width, WithinAbs(0.5 - 0.0867, 0.003));
        REQUIRE_THAT(top_inner_width, WithinAbs(0.45 - 0.0433, 0.003));
        REQUIRE_THAT(low_outer_reach, WithinAbs(9.75 + 0.078 + 0.026, 0.003));
        REQUIRE_THAT(low_inner_reach, WithinAbs(9.3 + 0.039, 0.003));
        REQUIRE_THAT(top_outer_reach, WithinAbs(9.75 - 0.0433 - 0.0433, 0.003));
        REQUIRE_THAT(top_inner_reach, WithinAbs(9.3 - 0.0217, 0.003));
        // The outer edge of the outer wall follows the outlines, the inner edge of the inner wall stays.
        REQUIRE_THAT(low_outer_reach + 0.5 * low_outer_width, WithinAbs(10.13, 0.005));
        REQUIRE_THAT(top_outer_reach + 0.5 * top_outer_width, WithinAbs(9.87, 0.005));
        REQUIRE_THAT(low_inner_reach - 0.5 * low_inner_width, WithinAbs(9.075, 0.005));
        REQUIRE_THAT(top_inner_reach - 0.5 * top_inner_width, WithinAbs(9.075, 0.005));
    }

    SECTION("an outer wall with no inner wall behind it takes what it can alone") {
        ExtrusionEntityCollection              perimeters = wall_loop(19.5);
        std::vector<ExtrusionEntityCollection> lower_passes;
        split_outer_wall_loops(perimeters, lower_passes, slices, slices.mid, Flow(0.5f, 0.1f, 0.4f), true, &inner_walls);
        REQUIRE(lower_passes.size() == 1);
        REQUIRE(lower_passes.front().entities.size() == 1);
        const auto [low_reach, low_width] = reach_and_width(lower_passes.front().entities[0]);
        REQUIRE_THAT(low_width, WithinAbs(0.6, 0.003));
        REQUIRE_THAT(low_reach, WithinAbs(9.8, 0.003));
    }

    SECTION("an inner wall behind a regular outer wall stays regular") {
        // The outer wall is not of the nominal width, so it is not printed in passes.
        ExtrusionEntityCollection perimeters = wall_loop(19.5);
        static_cast<ExtrusionLoop*>(perimeters.entities.front())->paths.front().width = 0.42f;
        perimeters.append(inner_wall_loop(18.6));
        std::vector<ExtrusionEntityCollection> lower_passes;
        split_outer_wall_loops(perimeters, lower_passes, slices, slices.mid, Flow(0.5f, 0.1f, 0.4f), true, &inner_walls);
        REQUIRE(lower_passes.size() == 1);
        REQUIRE(lower_passes.front().entities.empty());
    }
}

TEST_CASE("An inner wall around a hole shares the offset along its whole sides", "[OuterWallPasses]")
{
    // A square hole of 8 mm widening upwards by 0.13 mm a side per pass. The walls around it run along squares of
    // 8.5 mm and, at the distance the wall generator keeps, of 9.364 mm: the sides of the inner wall reach 0.68 mm
    // past the ends of the sides of the hole.
    OuterWallPassSlices slices;
    slices.mid       = { square(20., 8.) };
    slices.passes    = { { square(20., 7.74) }, { square(20., 8.26) } };
    slices.min_width = min_width;
    slices.max_width = max_width;
    slices.walls     = 2;
    const double                  first_spacing = 0.5 * (0.5 + 0.45) - 0.2 * (1. - 0.25 * PI);
    const OuterWallPassInnerWalls inner_walls { Flow(0.45f, 0.1f, 0.4f), first_spacing, 0.45 - 0.2 * (1. - 0.25 * PI) };
    ExtrusionEntityCollection     perimeters = wall_loop(8.5, true);
    perimeters.append(inner_wall_loop(8.5 + 2. * first_spacing, true));
    std::vector<ExtrusionEntityCollection> lower_passes;
    split_outer_wall_loops(perimeters, lower_passes, slices, slices.mid, Flow(0.5f, 0.1f, 0.4f), true, &inner_walls);
    REQUIRE(lower_passes.size() == 1);
    REQUIRE(lower_passes.front().entities.size() == 2);
    // The low pass of the inner wall takes its share of the 0.13 mm everywhere: nowhere it is left 0.45 mm wide.
    double narrowest = 1.;
    for (const ExtrusionPath &path : static_cast<const ExtrusionLoop*>(lower_passes.front().entities[1])->paths)
        narrowest = std::min(narrowest, double(path.width));
    REQUIRE(narrowest > 0.45 + 0.06);
}

TEST_CASE("Three walls around a hole share the offset along their whole sides", "[OuterWallPasses]")
{
    // A square hole of 8 mm widening upwards by 0.2 mm a side per pass: what three walls take inwards, 0.1 mm the
    // outer wall and 0.05 mm each inner wall. The sides of the third wall reach 1.09 mm past the ends of the sides
    // of the hole and 0.84 mm past those of the outer wall.
    OuterWallPassSlices slices;
    slices.mid       = { square(20., 8.) };
    slices.passes    = { { square(20., 7.6) }, { square(20., 8.4) } };
    slices.min_width = min_width;
    slices.max_width = max_width;
    slices.walls     = 3;
    const double                  spacing       = 0.45 - 0.2 * (1. - 0.25 * PI);
    const double                  first_spacing = 0.5 * (0.5 + 0.45) - 0.2 * (1. - 0.25 * PI);
    const OuterWallPassInnerWalls inner_walls { Flow(0.45f, 0.1f, 0.4f), first_spacing, spacing };
    ExtrusionEntityCollection     perimeters = wall_loop(8.5, true);
    perimeters.append(inner_wall_loop(8.5 + 2. * first_spacing, true));
    perimeters.append(inner_wall_loop(8.5 + 2. * (first_spacing + spacing), true, 2));
    std::vector<ExtrusionEntityCollection> lower_passes;
    split_outer_wall_loops(perimeters, lower_passes, slices, slices.mid, Flow(0.5f, 0.1f, 0.4f), true, &inner_walls);
    REQUIRE(lower_passes.size() == 1);
    REQUIRE(lower_passes.front().entities.size() == 3);
    // Outwards the walls share 0.2 mm as 0.1 to 0.15 to 0.15: the low pass of an inner wall gains 0.075 mm along
    // its sides and more at the corners, and is nowhere left 0.45 mm wide.
    for (const size_t wall : { size_t(1), size_t(2) }) {
        CAPTURE(wall);
        double narrowest = 1.;
        for (const ExtrusionPath &path : static_cast<const ExtrusionLoop*>(lower_passes.front().entities[wall])->paths)
            narrowest = std::min(narrowest, double(path.width));
        REQUIRE(narrowest > 0.45 + 0.06);
    }
    // The inner edge of the last wall stays: its middle moves by a half of what its width gains.
    const auto [low_reach, low_width] = reach_and_width(lower_passes.front().entities[2]);
    const auto [top_reach, top_width] = reach_and_width(perimeters.entities[2]);
    const double regular_reach = 0.5 * 8.5 + first_spacing + spacing;
    REQUIRE(low_width > 0.45);
    REQUIRE(top_width < 0.45 + 1e-6);
    REQUIRE_THAT(top_reach + 0.5 * top_width, WithinAbs(regular_reach + 0.225, 0.005));
}

TEST_CASE("Walls go up to the limits of the width together", "[OuterWallPasses]")
{
    // The outlines of the passes lie 0.35 mm off the middle one. Within 0.4 to 0.6 mm the two walls take 0.25 mm
    // outwards; the limits of 0.3 and 0.8 mm let the outer wall take 0.3 mm and the inner wall 0.35 mm.
    OuterWallPassSlices slices;
    slices.mid             = { square(20.) };
    slices.passes          = { { square(20.7) }, { square(19.3) } };
    slices.min_width       = min_width;
    slices.max_width       = max_width;
    slices.limit_min_width = 0.3;
    slices.limit_max_width = 0.8;
    slices.walls           = 2;
    const OuterWallPassInnerWalls inner_walls { Flow(0.45f, 0.1f, 0.4f), 0.45, 0.4 };
    ExtrusionEntityCollection     perimeters = wall_loop(19.5);
    perimeters.append(inner_wall_loop(18.6));
    std::vector<ExtrusionEntityCollection> lower_passes;
    split_outer_wall_loops(perimeters, lower_passes, slices, slices.mid, Flow(0.5f, 0.1f, 0.4f), true, &inner_walls);
    REQUIRE(lower_passes.size() == 1);
    REQUIRE(lower_passes.front().entities.size() == 2);
    // Of the 0.1 mm beyond the range the outer wall takes 0.2 / 0.4 and the inner wall the rest.
    const auto [low_outer_reach, low_outer_width] = reach_and_width(lower_passes.front().entities[0]);
    const auto [low_inner_reach, low_inner_width] = reach_and_width(lower_passes.front().entities[1]);
    REQUIRE_THAT(low_outer_width, WithinAbs(0.5 + 0.1 + 0.05, 0.003));
    REQUIRE_THAT(low_inner_width, WithinAbs(0.45 + 0.15 + 0.05, 0.003));
    REQUIRE_THAT(low_outer_reach + 0.5 * low_outer_width, WithinAbs(10.35, 0.005));
    REQUIRE_THAT(low_inner_reach - 0.5 * low_inner_width, WithinAbs(9.075, 0.005));
}

TEST_CASE("A pass may step back from the pass below it by the given step at most", "[OuterWallPasses]")
{
    // The upper pass steps back by 0.3 mm a side.
    const std::vector<ExPolygons> passes { { square(20.6) }, { square(20.) } };
    REQUIRE(outer_wall_pass_steps_fit(passes, nullptr, 0.35, true));
    REQUIRE_FALSE(outer_wall_pass_steps_fit(passes, nullptr, 0.25, true));
    // What lies over the top pass counts as well: an outline stepping back by 0.5 mm a side, or nothing at all.
    const ExPolygons above { square(19.) };
    const ExPolygons nothing;
    REQUIRE(outer_wall_pass_steps_fit(passes, &above, 0.55, true));
    REQUIRE_FALSE(outer_wall_pass_steps_fit(passes, &above, 0.35, true));
    REQUIRE_FALSE(outer_wall_pass_steps_fit(passes, &nothing, 0.35, true));
}

TEST_CASE("Holes count only when asked for", "[OuterWallPasses]")
{
    // The hole grows by 0.3 mm a side between the passes, the contour stays.
    const ExPolygons              mid { square(20., 8.) };
    const std::vector<ExPolygons> passes { { square(20., 7.4) }, { square(20., 8.6) } };
    REQUIRE_FALSE(outer_wall_passes_fit(mid, passes, min_width, max_width, true));
    REQUIRE(outer_wall_passes_fit(mid, passes, min_width, max_width, false));
}
