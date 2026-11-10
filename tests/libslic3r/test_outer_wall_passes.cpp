#include <catch2/catch_all.hpp>

#include <catch2/catch_test_macros.hpp>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/OuterWallPasses.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"

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

TEST_CASE("A part appearing or vanishing between the passes does not fit", "[OuterWallPasses]")
{
    const ExPolygons mid { square(20.) };
    ExPolygon island = square(2.);
    island.translate(Point(scaled<coord_t>(30.), 0));
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.), island }, { square(20.) } }, min_width, max_width, true));
    REQUIRE_FALSE(outer_wall_passes_fit(mid, { { square(20.) }, {} }, min_width, max_width, true));
}

TEST_CASE("Holes count only when asked for", "[OuterWallPasses]")
{
    // The hole grows by 0.3 mm a side between the passes, the contour stays.
    const ExPolygons              mid { square(20., 8.) };
    const std::vector<ExPolygons> passes { { square(20., 7.4) }, { square(20., 8.6) } };
    REQUIRE_FALSE(outer_wall_passes_fit(mid, passes, min_width, max_width, true));
    REQUIRE(outer_wall_passes_fit(mid, passes, min_width, max_width, false));
}
