#include <catch2/catch_all.hpp>
#include "libslic3r/GCode/PreciseSeam.hpp"

#include <algorithm>
#include <cmath>

using namespace Slic3r;

namespace {
Point mm(double x, double y) { return Point(scale_(x), scale_(y)); }

Polygon rectangle(double x0, double y0, double x1, double y1)
{
    // Modifier exteriors are CCW; holes are explicitly reversed in their fixtures.
    return Polygon(Points{mm(x0, y0), mm(x1, y0), mm(x1, y1), mm(x0, y1)});
}

void check_provenance(const Polygon &perimeter, const PreciseSeam::SegmentExtraction &result)
{
    REQUIRE(result.valid);
    for (const auto &segment : result.segments) {
        REQUIRE(segment.polyline.size() >= 2);
        REQUIRE(segment.edge_indices.size() + 1 == segment.polyline.size());
        CHECK(segment.length > 0.);
        const std::pair<PreciseSeam::PerimeterPosition, Point> endpoints[] = {
            {segment.begin, segment.polyline.points.front()}, {segment.end, segment.polyline.points.back()}};
        for (const auto &[position, point] : endpoints) {
            REQUIRE(position.edge_index < perimeter.size());
            CHECK(position.parameter >= 0.);
            CHECK(position.parameter < 1.);
            const Vec2d a = perimeter.points[position.edge_index].cast<double>();
            const Vec2d b = perimeter.points[(position.edge_index + 1) % perimeter.size()].cast<double>();
            const Vec2d reconstructed = a + position.parameter * (b - a);
            CHECK((point.cast<double>() - reconstructed).squaredNorm() <= 2.5);
        }
        for (size_t i = 0; i < segment.edge_indices.size(); ++i) {
            const size_t edge = segment.edge_indices[i];
            REQUIRE(edge < perimeter.size());
            const Vec2d a = perimeter.points[edge].cast<double>();
            const Vec2d direction = perimeter.points[(edge + 1) % perimeter.size()].cast<double>() - a;
            double previous = -1.;
            for (const Point &point : {segment.polyline.points[i], segment.polyline.points[i + 1]}) {
                const Vec2d offset = point.cast<double>() - a;
                const double t = offset.dot(direction) / direction.squaredNorm();
                CHECK(t >= 0.);
                CHECK(t <= 1.);
                CHECK(t >= previous);
                CHECK((offset - t * direction).squaredNorm() <= 2.5);
                previous = t;
            }
        }
    }
}
} // namespace

TEST_CASE("A crossing modifier extracts both perimeter intervals without a body chord", "[PreciseSeam][SegmentExtraction]")
{
    const bool reverse = GENERATE(false, true);
    Polygon perimeter = rectangle(0, 0, 20, 20);
    if (reverse) perimeter.reverse();
    const Points original = perimeter.points;
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(8, -2, 12, 22))});
    check_provenance(perimeter, result);
    CHECK(perimeter.points == original);
    REQUIRE(result.segments.size() == 2);
    CHECK_FALSE(result.full_containment);
    std::vector<coord_t> sides;
    for (const auto &segment : result.segments) {
        CHECK_THAT(unscale<double>(segment.length), Catch::Matchers::WithinAbs(4., 1e-6));
        CHECK(segment.polyline.points.front().y() == segment.polyline.points.back().y());
        sides.push_back(segment.polyline.points.front().y());
    }
    std::sort(sides.begin(), sides.end());
    CHECK(sides == std::vector<coord_t>{mm(0, 0).y(), mm(0, 20).y()});
}

TEST_CASE("A corner interval remains connected when the contour origin changes", "[PreciseSeam][SegmentExtraction]")
{
    const size_t origin = GENERATE(size_t(0), size_t(1), size_t(2), size_t(3));
    Polygon perimeter = rectangle(0, 0, 20, 20);
    std::rotate(perimeter.points.begin(), perimeter.points.begin() + origin, perimeter.points.end());
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(-2, -2, 4, 4))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 1);
    const auto &segment = result.segments.front();
    CHECK(segment.polyline.points.front() == mm(0, 4));
    CHECK(segment.polyline.points.back() == mm(4, 0));
    CHECK_THAT(unscale<double>(segment.length), Catch::Matchers::WithinAbs(8., 1e-6));
    CHECK_FALSE(result.full_containment);
}

TEST_CASE("Collinear perimeter vertices retain their original edge provenance", "[PreciseSeam][SegmentExtraction]")
{
    const Polygon perimeter(Points{mm(0, 0), mm(2, 0), mm(4, 0), mm(8, 0), mm(20, 0), mm(20, 20), mm(0, 20)});
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(1, -2, 13, 2))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 1);
    const auto &segment = result.segments.front();
    CHECK(segment.polyline.points.front() == mm(1, 0));
    CHECK(segment.polyline.points.back() == mm(13, 0));
    CHECK(segment.edge_indices == std::vector<size_t>{0, 1, 2, 3});
    CHECK_THAT(unscale<double>(segment.length), Catch::Matchers::WithinAbs(12., 1e-6));
}

TEST_CASE("Interior vertices retain their sequence between two cut endpoints", "[PreciseSeam][SegmentExtraction]")
{
    const bool reverse = GENERATE(false, true);
    const size_t origin = GENERATE(size_t(0), size_t(5));
    Polygon perimeter;
    // A zigzag makes skipped or misbound interior edges observable in the arc length.
    for (int x = 0; x <= 10; ++x)
        perimeter.points.push_back(mm(x, x % 2));
    perimeter.points.push_back(mm(10, 10));
    perimeter.points.push_back(mm(0, 10));
    if (reverse) perimeter.reverse();
    std::rotate(perimeter.points.begin(), perimeter.points.begin() + origin, perimeter.points.end());
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(2.5, -2, 8.5, 3))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 1);
    const auto &segment = result.segments.front();
    Points expected{mm(2.5, .5)};
    for (int x = 3; x <= 8; ++x)
        expected.push_back(mm(x, x % 2));
    expected.push_back(mm(8.5, .5));
    if (reverse) std::reverse(expected.begin(), expected.end());
    CHECK(segment.polyline.points == expected);
    CHECK_THAT(unscale<double>(segment.length), Catch::Matchers::WithinAbs(6. * std::sqrt(2.), 1e-6));
}

TEST_CASE("Neighboring vertices distinguish repeated anchors on different lobes", "[PreciseSeam][SegmentExtraction]")
{
    const bool reverse = GENERATE(false, true);
    // Both lobes visit the origin, but their adjacent edges lead to different vertices.
    Polygon perimeter(Points{mm(0, 0), mm(4, 0), mm(4, 4), mm(0, 4),
                             mm(0, 0), mm(-4, 0), mm(-4, -4), mm(0, -4)});
    if (reverse) perimeter.reverse();
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(-1, -5, 1, 5))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 2);
    std::vector<size_t> edges;
    for (const auto &segment : result.segments) {
        CHECK_THAT(unscale<double>(segment.length), Catch::Matchers::WithinAbs(6., 1e-6));
        edges.insert(edges.end(), segment.edge_indices.begin(), segment.edge_indices.end());
    }
    std::sort(edges.begin(), edges.end());
    CHECK(edges == std::vector<size_t>{0, 2, 3, 4, 6, 7});
}

TEST_CASE("Modifier holes subtract coverage while separate components add intervals", "[PreciseSeam][SegmentExtraction]")
{
    const Polygon perimeter = rectangle(0, 0, 20, 20);
    ExPolygon area(rectangle(1, -3, 10, 3));
    area.holes.push_back(rectangle(4, -1, 7, 1));
    area.holes.back().reverse();
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter,
        {area, ExPolygon(rectangle(14, -3, 18, 3))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 3);
    const std::vector<Point> starts{mm(1, 0), mm(7, 0), mm(14, 0)};
    const std::vector<Point> ends{mm(4, 0), mm(10, 0), mm(18, 0)};
    for (size_t i = 0; i < 3; ++i) {
        CHECK(result.segments[i].polyline.points.front() == starts[i]);
        CHECK(result.segments[i].polyline.points.back() == ends[i]);
    }
    CHECK_FALSE(result.full_containment);
}

TEST_CASE("Full coverage is distinct from an empty or point-only intersection", "[PreciseSeam][SegmentExtraction]")
{
    const Polygon perimeter = rectangle(0, 0, 20, 20);
    const int scenario = GENERATE(0, 1, 2, 3, 4);
    ExPolygons modifier;
    if (scenario == 0) modifier = {ExPolygon(rectangle(-2, -2, 22, 22))};
    if (scenario == 1) modifier = {ExPolygon(rectangle(2, 2, 4, 4))};
    if (scenario == 2) modifier = {ExPolygon(rectangle(30, 30, 40, 40))};
    if (scenario == 3) modifier = {ExPolygon(rectangle(20, 20, 25, 25))};
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, modifier);
    check_provenance(perimeter, result);
    CHECK(result.full_containment == (scenario == 0));
    if (scenario == 0) {
        REQUIRE(result.segments.size() == 1);
        CHECK_THAT(unscale<double>(result.segments.front().length), Catch::Matchers::WithinAbs(80., 1e-6));
        CHECK(result.segments.front().edge_indices.size() == 4);
    } else
        CHECK(result.segments.empty());
}

TEST_CASE("Repeated visits to a coordinate stay on their original perimeter edges", "[PreciseSeam][SegmentExtraction]")
{
    // Two visits to the origin belong to different lobes, not to one shared vertex.
    const Polygon perimeter(Points{mm(0, 0), mm(4, 0), mm(4, 4), mm(0, 0), mm(-4, 0), mm(-4, -4)});
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(-1, -1, 1, 1))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 2);
    std::vector<size_t> edges;
    for (const auto &segment : result.segments)
        edges.insert(edges.end(), segment.edge_indices.begin(), segment.edge_indices.end());
    std::sort(edges.begin(), edges.end());
    CHECK(edges == std::vector<size_t>{0, 2, 3, 5});
}

TEST_CASE("Two-point fragments accept the first matching perimeter edge", "[PreciseSeam][SegmentExtraction]")
{
    // Policy: do not search for duplicate bindings on overlapping source edges.
    const Polygon perimeter(Points{mm(0, 0), mm(10, 0), mm(0, 0), mm(0, 10)});
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(2, -1, 8, 1))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 1);
    CHECK(result.segments[0].edge_indices == std::vector<size_t>{0});
    CHECK(result.segments[0].polyline.points.front() == mm(2, 0));
    CHECK(result.segments[0].polyline.points.back() == mm(8, 0));
    CHECK_THAT(unscale<double>(result.segments[0].length), Catch::Matchers::WithinAbs(6., 1e-6));
}

TEST_CASE("Unnormalized perimeter input is reported instead of silently losing coverage", "[PreciseSeam][SegmentExtraction]")
{
    const int scenario = GENERATE(0, 1, 2);
    Polygon perimeter;
    if (scenario == 1) perimeter.points = {mm(0, 0), mm(20, 0)};
    if (scenario == 2) perimeter.points = {mm(0, 0), mm(20, 0), mm(20, 0), mm(0, 20)};
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(-2, -2, 22, 22))});
    CHECK_FALSE(result.valid);
    CHECK(result.segments.empty());
    CHECK_FALSE(result.full_containment);
}

TEST_CASE("Segment lengths measure diagonal arcs rather than squared distances", "[PreciseSeam][SegmentExtraction]")
{
    const Polygon perimeter(Points{mm(0, 0), mm(10, 10), mm(0, 10)});
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(2, -1, 6, 11))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 2);
    CHECK(result.segments[0].polyline.points.front() == mm(2, 2));
    CHECK(result.segments[0].polyline.points.back() == mm(6, 6));
    CHECK_THAT(unscale<double>(result.segments[0].length), Catch::Matchers::WithinAbs(std::sqrt(32.), 1e-6));
    CHECK_THAT(unscale<double>(result.segments[1].length), Catch::Matchers::WithinAbs(4., 1e-6));
}

TEST_CASE("An endpoint at an original vertex uses its outgoing edge including vertex zero", "[PreciseSeam][SegmentExtraction]")
{
    const size_t origin = GENERATE(size_t(0), size_t(1));
    Polygon perimeter(Points{mm(0, 0), mm(5, 5), mm(0, 10)});
    std::rotate(perimeter.points.begin(), perimeter.points.begin() + origin, perimeter.points.end());
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {ExPolygon(rectangle(-1, -1, 6, 5))});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 1);
    const auto &segment = result.segments.front();
    CHECK(segment.polyline.points.back() == mm(5, 5));
    CHECK(segment.end.edge_index == 1 - origin);
    CHECK_THAT(segment.end.parameter, Catch::Matchers::WithinAbs(0., 1e-12));
}

TEST_CASE("A one-unit uncovered gap is not bridged by the projection tolerance", "[PreciseSeam][SegmentExtraction]")
{
    // These are scaled integer units, deliberately smaller than the projection tolerance.
    const Polygon perimeter(Points{Point(0, 0), Point(20, 0), Point(20, 20), Point(0, 20)});
    const ExPolygons modifier{
        ExPolygon(Polygon(Points{Point(3, -2), Point(8, -2), Point(8, 2), Point(3, 2)})),
        ExPolygon(Polygon(Points{Point(9, -2), Point(14, -2), Point(14, 2), Point(9, 2)}))};
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, modifier);
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 2);
    CHECK(result.segments[0].polyline.points.back() == Point(8, 0));
    CHECK(result.segments[1].polyline.points.front() == Point(9, 0));
}

TEST_CASE("Distant modifier areas do not change nearby coverage or detach its holes", "[PreciseSeam][SegmentExtraction]")
{
    const bool reverse_areas = GENERATE(false, true);
    const Polygon perimeter = rectangle(0, 0, 20, 20);
    ExPolygon nearby(rectangle(1, -3, 10, 3));
    nearby.holes.push_back(rectangle(4, -1, 7, 1));
    nearby.holes.back().reverse();
    // A rejected area keeps its own hole; neither may affect the nearby area.
    ExPolygon distant(rectangle(100, 100, 120, 120));
    distant.holes.push_back(rectangle(105, 105, 115, 115));
    distant.holes.back().reverse();
    ExPolygons modifier{distant, nearby, ExPolygon(rectangle(14, -3, 18, 3))};
    if (reverse_areas) std::reverse(modifier.begin(), modifier.end());
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, modifier);
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 3);
    const Points starts{mm(1, 0), mm(7, 0), mm(14, 0)};
    const Points ends{mm(4, 0), mm(10, 0), mm(18, 0)};
    for (size_t i = 0; i < starts.size(); ++i) {
        CHECK(result.segments[i].polyline.points.front() == starts[i]);
        CHECK(result.segments[i].polyline.points.back() == ends[i]);
        CHECK(result.segments[i].edge_indices == std::vector<size_t>{0});
    }
    CHECK_FALSE(result.full_containment);
}

TEST_CASE("Rounded intersections on an inclined edge retain their original edge", "[PreciseSeam][SegmentExtraction]")
{
    const bool reverse = GENERATE(false, true);
    // Grid intersections (3, 0.9) and (7, 2.1) round to (3, 1) and (7, 2).
    Polygon perimeter(Points{Point(0, 0), Point(10, 3), Point(10, 20), Point(0, 20)});
    if (reverse) perimeter.reverse();
    const ExPolygon area(Polygon(Points{Point(3, -2), Point(7, -2), Point(7, 5), Point(3, 5)}));
    const auto result = PreciseSeam::extract_perimeter_segments(perimeter, {area});
    check_provenance(perimeter, result);
    REQUIRE(result.segments.size() == 1);
    const auto &segment = result.segments.front();
    CHECK(segment.polyline.points.front() == (reverse ? Point(7, 2) : Point(3, 1)));
    CHECK(segment.polyline.points.back() == (reverse ? Point(3, 1) : Point(7, 2)));
    CHECK(segment.edge_indices == std::vector<size_t>{reverse ? size_t(2) : size_t(0)});
    CHECK_THAT(segment.length, Catch::Matchers::WithinAbs(std::sqrt(17.), 1e-9));
    CHECK_FALSE(result.full_containment);
}
