#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <string>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ArcFitter.hpp"
#include "libslic3r/Circle.hpp"
#include "libslic3r/GCode/AISeamPlacement.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r/libslic3r.h"
#include "test_utils.hpp"

using namespace Slic3r;
using namespace Slic3r::AISeamPlacement;

namespace {
ExtrusionLoop square()
{
    // Non-planar first edge distinguishes interpolation from copying either endpoint's Z.
    ExtrusionPath path(erExternalPerimeter, 0.08, 0.4f, 0.2f);
    path.z_contoured = true;
    path.polyline.points = {Point3(0, 0, 0), Point3(scale_(10.), 0., scale_(1.)),
        Point3(scale_(10.), scale_(10.), 0.), Point3(0., scale_(10.), 0.), Point3(0, 0, 0)};
    return ExtrusionLoop(path);
}

Json response()
{
    return {{"format", "orca-ai-seams-response"}, {"version", 1}, {"export_id", "example"},
        {"export_digest", std::string(64, 'a')},
        {"seams", {{{"occurrence_id", "loop/i0"}, {"x_mm", 0.}, {"y_mm", 0.}, {"edge", {{0., 0.}, {10., 0.}}}}}}};
}
} // namespace

TEST_CASE("Approved seams retain submicron positions and interpolate contour Z", "[AISeamPlacement]")
{
    ExtrusionLoop loop = square();
    const double x = 2. * SCALING_FACTOR; // Inside the standard 0.0015 mm vertex snap radius.
    const Projection projected = project(loop, Vec2d(x, 0.));
    split_exact(loop, projected);
    REQUIRE(loop.first_point() == Point(2, 0));
    REQUIRE(loop.paths.front().first_point3().z() == coord_t(std::llround(0.2)));
    REQUIRE(loop.paths.front().z_contoured);
    REQUIRE(loop.paths.front().role() == erExternalPerimeter);
    REQUIRE_THAT(loop.paths.front().mm3_per_mm, Catch::Matchers::WithinAbs(0.08, 1e-12));

    loop = square();
    split_exact(loop, project(loop, Vec2d(2.5, 0.)));
    REQUIRE(loop.paths.front().first_point3().z() == scale_(0.25));
}

TEST_CASE("Approved overhang seams do not move to supported paths", "[AISeamPlacement]")
{
    ExtrusionLoop loop = square();
    loop.paths.front().set_extrusion_role(erOverhangPerimeter);
    loop.inset_idx = 0; // Fully overhanging outer walls are still exportable.
    REQUIRE(is_external(loop));
    split_exact(loop, project(loop, Vec2d(5., 0.)));
    REQUIRE(loop.first_point() == Point(scale_(5.), 0.));
    REQUIRE(loop.paths.front().role() == erOverhangPerimeter);
    loop.inset_idx = 1;
    REQUIRE_FALSE(is_external(loop));
}

TEST_CASE("Projection reports distances without snapping nearby vertices", "[AISeamPlacement]")
{
    const ExtrusionLoop loop = square();
    const Projection projection = project(loop, Vec2d(5., 0.01));
    REQUIRE_THAT(projection.distance_mm, Catch::Matchers::WithinAbs(0.01, 1e-12));
    REQUIRE(projection.point == Point(scale_(5.), 0.));
    REQUIRE_THAT(projection.s_mm, Catch::Matchers::WithinAbs(5., 1e-12));
    REQUIRE_THROWS(project(loop, Vec2d(std::numeric_limits<double>::infinity(), 0.)));
}

TEST_CASE("A short response edge cannot open a chain around the loop", "[AISeamPlacement]")
{
    // Review counterexample: A and B are within the old 0.05 mm vertex radius of each other.
    ExtrusionPath path(erExternalPerimeter, 0.08, 0.4f, 0.2f);
    path.polyline.points = {Point3(0., 0., 0.), Point3(scale_(0.01), scale_(0.01), 0.), Point3(scale_(10.), 0., 0.),
        Point3(scale_(10.), scale_(10.), 0.), Point3(0., scale_(10.), 0.), Point3(0., 0., 0.)};
    const ExtrusionLoop loop(path);
    const Vec2d a(0., 0.), b(0.01, 0.01);
    const auto far = project_on_edge(loop, Vec2d(5., 10.), a, b);
    REQUIRE(far.has_value());
    REQUIRE(far->distance_mm > 1.); // Far from A-B, never moved onto D-E.
    const auto middle = project_on_edge(loop, Vec2d(0.005, 0.005), a, b);
    REQUIRE(middle.has_value());
    REQUIRE(middle->distance_mm < 1e-6);
    const auto swapped = project_on_edge(loop, Vec2d(0.005, 0.005), b, a);
    REQUIRE(swapped.has_value());
    REQUIRE(swapped->distance_mm < 1e-6);
    const auto closing = project_on_edge(loop, Vec2d(0., 5.), Vec2d(0., 10.), a); // E-A closes the loop.
    REQUIRE(closing.has_value());
    REQUIRE(closing->distance_mm < 1e-6);
    REQUIRE(closing->segment_idx == 4);
    REQUIRE_FALSE(project_on_edge(loop, Vec2d(5., 10.), a, Vec2d(10., 10.)).has_value()); // A-D is not a segment.
}

TEST_CASE("A fitting swapped edge beats a matching but unfitting path-order edge", "[AISeamPlacement]")
{
    // Review round 3: rounded ends of the thin loop match A-B in path order and C-D swapped.
    ExtrusionPath path(erExternalPerimeter, 0.08, 0.4f, 0.2f);
    path.polyline.points = {Point3(0., scale_(0.004), 0.), Point3(scale_(10.), scale_(0.004), 0.),
        Point3(scale_(10.), scale_(0.014), 0.), Point3(0., scale_(0.014), 0.), Point3(0., scale_(0.004), 0.)};
    const ExtrusionLoop loop(path);
    const auto projection = project_on_edge(loop, Vec2d(5., 0.06), Vec2d(0., 0.01), Vec2d(10., 0.01));
    REQUIRE(projection.has_value());
    REQUIRE(projection->segment_idx == 2); // C-D
    REQUIRE(projection->distance_mm < 0.05);
}

TEST_CASE("Responses reject duplicate IDs and nonfinite coordinates", "[AISeamPlacement]")
{
    Json value = response();
    REQUIRE_NOTHROW(validate_response(value));
    value["seams"].push_back(value["seams"][0]);
    REQUIRE_THROWS(validate_response(value));
    value = response();
    value["seams"][0]["x_mm"] = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS(validate_response(value));
    value = response();
    value["seams"][0]["y_mm"] = "0.0";
    REQUIRE_THROWS(validate_response(value));
    value = response();
    value["seams"][0].erase("edge");
    REQUIRE_THROWS(validate_response(value));
    value = response();
    value["seams"][0]["edge"][1] = {10.};
    REQUIRE_THROWS(validate_response(value));
}

TEST_CASE("Response parsing rejects duplicate JSON object members", "[AISeamPlacement]")
{
    ScopedTemporaryFile file(".json");
    // Repeated keys must not be silently resolved by last-member-wins parsing.
    const std::string text = "{\"format\":\"orca-ai-seams-response\",\"format\":\"other\"}";
    REQUIRE_FALSE(bool(write_file_atomically(file.string(), text)));
    REQUIRE_THROWS(read_response(file.string()));
}

TEST_CASE("Exact splitting retains path attributes and linear fitting metadata", "[AISeamPlacement]")
{
    ExtrusionLoop loop = square();
    auto &path = loop.paths.front();
    path.overhang_degree = 0.25;
    path.curve_degree = 7;
    path.smooth_speed = 42.;
    path.polyline.fitting_result.push_back(
        {0, path.polyline.points.size() - 1, EMovePathType::Linear_move, ArcSegment()});
    split_exact(loop, project(loop, Vec2d(2.5, 0.)));
    bool has_fitting = false;
    for (const auto &part : loop.paths) {
        REQUIRE(part.z_contoured);
        REQUIRE(part.curve_degree == 7);
        REQUIRE_THAT(part.overhang_degree, Catch::Matchers::WithinAbs(0.25, 1e-12));
        REQUIRE_THAT(part.smooth_speed, Catch::Matchers::WithinAbs(42., 1e-12));
        REQUIRE_THAT(part.width, Catch::Matchers::WithinAbs(0.4, 1e-6));
        REQUIRE_THAT(part.height, Catch::Matchers::WithinAbs(0.2, 1e-6));
        for (const auto &fit : part.polyline.fitting_result) {
            has_fitting = true;
            REQUIRE(fit.start_point_index <= fit.end_point_index);
            REQUIRE(fit.end_point_index < part.polyline.points.size());
        }
    }
    REQUIRE(has_fitting);
}

TEST_CASE("Exact splitting retains the selected visit at repeated coordinates", "[AISeamPlacement]")
{
    ExtrusionLoop loop = square();
    loop.paths.front().polyline.points = {Point3(0, 0, 0), Point3(10, 0, 0), Point3(0, 0, 0),
        Point3(0, 10, 0), Point3(0, 0, 0)};
    Projection selected;
    selected.segment_idx = 2;
    selected.point = Point(0, 0);
    split_exact(loop, selected);
    REQUIRE(loop.paths.front().polyline.points[1].to_point() == Point(0, 10));
}
