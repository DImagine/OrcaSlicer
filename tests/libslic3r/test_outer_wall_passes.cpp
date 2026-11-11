#include <catch2/catch_all.hpp>

#include <catch2/catch_test_macros.hpp>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/OuterWallPasses.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/libslic3r.h"

#include <cstddef>
#include <utility>
#include <vector>

using namespace Slic3r;

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

// Outer wall of the nominal width as one loop along a square of the given side, counter-clockwise.
ExtrusionEntityCollection wall_loop(double side)
{
    const coord_t half = scaled<coord_t>(0.5 * side);
    ExtrusionPath path(erExternalPerimeter, 0.1, float(0.5 * (min_width + max_width)), 0.2f);
    path.polyline = Polyline3(Polyline(Points{ Point(-half, -half), Point(half, -half), Point(half, half), Point(-half, half),
                                               Point(-half, -half) }));
    ExtrusionLoop loop(path);
    loop.inset_idx = 0;
    ExtrusionEntityCollection out;
    out.append(loop);
    return out;
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
    // Number of the low passes made of a wall loop running along a square of 19.5 mm, for the outlines of the middle
    // of the layer and of its two passes.
    auto low_passes = [](const ExPolygon &mid, const ExPolygon &low, const ExPolygon &high) {
        OuterWallPassSlices slices;
        slices.mid       = { mid };
        slices.passes    = { { low }, { high } };
        slices.min_width = min_width;
        slices.max_width = max_width;
        ExtrusionEntityCollection              perimeters = wall_loop(19.5);
        std::vector<ExtrusionEntityCollection> lower_passes;
        split_outer_wall_loops(perimeters, lower_passes, slices, Flow(0.5f, 0.1f, 0.4f), true);
        return lower_passes.empty() ? size_t(0) : lower_passes.front().entities.size();
    };
    // A solid square of 20 mm narrowing upwards: the loop runs half its width inside the outline.
    REQUIRE(low_passes(square(20.), square(20.05), square(19.95)) == 1);
    // A frame 0.2 mm thick, the wall thickened inwards by an XY compensation of its hole the outlines do not
    // carry. The loop runs in the hole of the outlines, and their nearest edge is not the one it follows.
    REQUIRE(low_passes(square(20., 19.6), square(20.05, 19.65), square(19.95, 19.55)) == 0);
}

TEST_CASE("Holes count only when asked for", "[OuterWallPasses]")
{
    // The hole grows by 0.3 mm a side between the passes, the contour stays.
    const ExPolygons              mid { square(20., 8.) };
    const std::vector<ExPolygons> passes { { square(20., 7.4) }, { square(20., 8.6) } };
    REQUIRE_FALSE(outer_wall_passes_fit(mid, passes, min_width, max_width, true));
    REQUIRE(outer_wall_passes_fit(mid, passes, min_width, max_width, false));
}
