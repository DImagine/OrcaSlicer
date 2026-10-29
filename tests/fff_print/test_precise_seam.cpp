#include <catch2/catch_all.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "test_helpers.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/GCode/PreciseSeam.hpp"
#include "libslic3r/GCode/PreciseSeamRules.hpp"
#include "libslic3r/GCode/PreciseSeamRulesApply.hpp"
#include "libslic3r/GCode/PreciseSeamRulesExport.hpp"
#include "libslic3r/GCode/ToolOrdering.hpp"
#include "libslic3r/GCode/PreciseSeamRulesGeometry.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"

#include <algorithm>
#include <functional>
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/PrintConfig.hpp"
#include <vector>
#include <optional>
#include <utility>
#include "libslic3r/libslic3r.h"
#include "libslic3r/GCode/SeamPlacer.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/ExPolygon.hpp"
#include <cstddef>
#include <string>
#include <string_view>
#include <cstdlib>

using namespace Slic3r;

namespace {
Point mm(double x, double y) { return Point(scale_(x), scale_(y)); }

Polygon rectangle(double x0, double y0, double x1, double y1)
{
    // Counterclockwise contours match the modifier-slice cache contract.
    return Polygon(Points{mm(x0, y0), mm(x1, y0), mm(x1, y1), mm(x0, y1)});
}

struct SeamFixture {
    Model model;
    Print print;
    Model modifiers;
    Layer *layer = nullptr;
    PreciseSeam::ModifierRegionsCache cache;

    explicit SeamFixture(int raft_layers = 0)
    {
        // Only the layer/PrintObject context is needed; clipping uses explicit cached slices below.
        Test::init_print({Test::cube(20)}, print, model, {{"raft_layers", std::to_string(raft_layers)}});
        REQUIRE(print.objects().size() == 1);
        PrintObject *object = print.get_object(0);
        const auto &slicing = object->slicing_parameters();
        // IDs and print heights include the raft; mesh slicing heights remain object-relative.
        layer = object->add_layer(int(slicing.raft_layers()), 0.2, slicing.object_print_z_min + 0.2, 0.1);
        modifiers.add_object();
    }

    const ModelVolume *add(ModelVolumeType type, Polygons slices)
    {
        // These volumes own cache keys; mesh slicing is deliberately outside this geometry fixture.
        ModelVolume *volume = modifiers.objects.front()->add_volume(Test::cube(1));
        volume->set_type(type);
        ExPolygons regions;
        for (Polygon &slice : slices)
            regions.emplace_back(std::move(slice));
        cache.emplace(volume, PreciseSeam::prepare_modifier_slices({std::move(regions)}));
        return volume;
    }

    const ModelVolume *add_regions(ModelVolumeType type, ExPolygons regions)
    {
        // Supply structured slices directly, without reconstructing holes from flat contours.
        const ModelVolume *volume = add(type, {});
        REQUIRE(volume->is_precise_seam());
        cache.at(volume) = PreciseSeam::prepare_modifier_slices({std::move(regions)});
        return volume;
    }
};

void check_square_boundary(const Polygon &polygon)
{
    // Wrong insertion edges can retrace a side without changing area: check length as well.
    CHECK_THAT(unscale<double>(polygon.length()), Catch::Matchers::WithinAbs(80.0, 0.00001));
    for (const Point &p : polygon.points) {
        CAPTURE(p.x(), p.y());
        CHECK(p.x() >= scale_(0));
        CHECK(p.x() <= scale_(20));
        CHECK(p.y() >= scale_(0));
        CHECK(p.y() <= scale_(20));
        const bool on_boundary = p.x() == 0 || p.x() == scale_(20) || p.y() == 0 || p.y() == scale_(20);
        CHECK(on_boundary);
    }
}

void require_vertex(const Polygon &polygon, const Point &point)
{
    // Integer coordinates make the micron transition helpers exact on these axis-aligned edges.
    CAPTURE(point.x(), point.y());
    REQUIRE(std::find(polygon.points.begin(), polygon.points.end(), point) != polygon.points.end());
}

std::vector<SeamPlacerImpl::EnforcedBlockedSeamPoint> weak_candidate_types(
    const Polygon &polygon, const std::vector<PreciseSeam::WeakModifierSegment> &segments,
    const std::function<SeamPlacerImpl::EnforcedBlockedSeamPoint(const Point &)> &painted = {},
    PreciseSeam::PreciseSeamWarnings *warnings = nullptr)
{
    // Use the same coordinate conversion as production when applying prepared boundaries.
    // Optional painted types stand in for seam painting, which production assigns before weak zones.
    PrintObjectSeamData::LayerSeams candidates;
    candidates.perimeters.emplace_back();
    auto &loop = candidates.perimeters.back();
    loop.start_index = 0;
    for (const Point &point : polygon.points) {
        const Vec2f position = unscale(point).cast<float>();
        candidates.points.emplace_back(Vec3f(position.x(), position.y(), 0), loop, 0,
                                       painted ? painted(point) : SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral);
    }
    loop.end_index = candidates.points.size();
    bool enforced = false;
    PreciseSeam::apply_weak_modifiers_to_perimeter(segments, candidates, loop, enforced, warnings);
    std::vector<SeamPlacerImpl::EnforcedBlockedSeamPoint> types;
    for (const auto &candidate : candidates.points)
        types.push_back(candidate.type);
    return types;
}
} // namespace

TEST_CASE("A simple clipped interval has consistent geometry before and after weak boundary insertion", "[PreciseSeam][SegmentExtraction]")
{
    const bool corner = GENERATE(false, true);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Polygon cut = corner ? rectangle(-2, -2, 4, 4) : rectangle(2, -2, 8, 2);
    const Point expected_begin = corner ? mm(0, 4) : mm(2, 0);
    const Point expected_end = corner ? mm(4, 0) : mm(8, 0);
    const auto extracted = PreciseSeam::extract_perimeter_segments(
        PreciseSeam::PreparedPerimeter(perimeter), PreciseSeam::prepare_modifier_regions(ExPolygons{ExPolygon(cut)}), ModelVolumeType::PRECISE_SEAM_ENFORCED);
    REQUIRE(extracted.valid);
    REQUIRE(extracted.segments.size() == 1);
    CHECK(extracted.segments.front().polyline.points.front() == expected_begin);
    CHECK(extracted.segments.front().polyline.points.back() == expected_end);
    // Boundary insertion must preserve the extractor endpoints on this analytic fixture.
    const auto *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {cut});
    const auto applied = PreciseSeam::collect_weak_modifier_segments({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(applied.size() == 1);
    CHECK(applied.front().left_point == expected_begin);
    CHECK(applied.front().right_point == expected_end);
    check_square_boundary(perimeter);
}

TEST_CASE("Unsuccessful strong modifiers share perimeter preparation with weak processing", "[PreciseSeam][Regression]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Points original = perimeter.points;
    const auto *empty = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {});
    const auto *distant = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(30, 30, 40, 40)});
    const auto *covering = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(-2, -2, 22, 22)});
    const auto *blocked = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(2, -2, 8, 2)});
    const auto *neutral = fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, {rectangle(4, -2, 6, 2)});
    const PreciseSeam::PreparedPerimeter prepared(perimeter);
    REQUIRE(prepared.valid);
    CHECK(prepared.bounds.min == mm(0, 0));
    CHECK(prepared.bounds.max == mm(20, 20));
    REQUIRE(prepared.line.size() == original.size() + 1);
    CHECK(prepared.line.points.front() == prepared.line.points.back());
    const Point *line_storage = prepared.line.points.data();
    PreciseSeam::PreciseSeamWarnings warnings;
    CHECK_FALSE(PreciseSeam::insert_strong_seam_point(
        {empty, distant, covering}, perimeter, prepared, fixture.layer, fixture.cache, &warnings).has_value());
    CHECK(perimeter.points == original);
    CHECK(prepared.line.points.data() == line_storage);
    CHECK(warnings.full_containment.load());

    // Weak receives the same preparation; neither consumer needs to reconstruct the clipping line.
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {blocked, neutral}, perimeter, prepared, fixture.layer, fixture.cache, &warnings);
    REQUIRE(segments.size() == 2);
    // Insertion has now changed the polygon: do not use prepared for any further extraction.
    const auto types = weak_candidate_types(perimeter, segments);
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &p = perimeter[i];
        const bool in_outer = p.y() == 0 && p.x() >= scale_(2.) && p.x() <= scale_(8.);
        const bool in_inner = p.x() >= scale_(4.) && p.x() <= scale_(6.);
        const auto expected = in_outer && !in_inner ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                                                     SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral;
        CHECK(types[i] == expected);
    }
    CHECK(warnings.failed_fragments.load() == 0);
    check_square_boundary(perimeter);
}

TEST_CASE("Structured modifier slices keep holes with their component and preserve sliced area", "[PreciseSeam][SegmentExtraction]")
{
    SeamFixture fixture;
    TriangleMesh shell = Test::cube(10);
    TriangleMesh cavity = Test::cube(6);
    cavity.translate(2., 2., -1.);
    cavity.flip_triangles();
    shell.merge(cavity);
    TriangleMesh island = Test::cube(2);
    island.translate(30., 0., 0.);
    shell.merge(island);
    // The layer intersects an annulus and a separate island, with areas 100-36 and 4 mm^2.
    ModelVolume *volume = fixture.modifiers.objects.front()->add_volume(shell);
    const bool mirrored = GENERATE(false, true);
    // A negative determinant must preserve exterior/hole winding and the sliced area.
    if (mirrored)
        volume->set_mirror(Vec3d(-1., 1., 1.));
    PrintObject *object = fixture.print.get_object(0);
    PreciseSeam::ModifierRegionsCache cache;
    cache.emplace(volume, PreciseSeam::prepare_modifier_slices(object->slice_single_volume_regions(volume)));
    const auto &layers = cache.at(volume);
    REQUIRE(layers.size() == 1);
    REQUIRE(layers.front().size() == 2);
    size_t holes = 0;
    double area = 0.;
    for (const auto &cached_region : layers.front()) {
        const ExPolygon &region = cached_region.polygon;
        holes += region.holes.size();
        CHECK(region.contour.is_counter_clockwise());
        for (const Polygon &hole : region.holes)
            CHECK(hole.is_clockwise());
        area += region.area();
    }
    CHECK(holes == 1);
    CHECK_THAT(area / double(scale_(1.)) / double(scale_(1.)), Catch::Matchers::WithinAbs(68., 1e-4));
}

TEST_CASE("Modifier slices above a raft use object layer indices for strong and weak seams", "[PreciseSeam][Regression]")
{
    const auto mode = GENERATE(ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_CENTER,
                              ModelVolumeType::PRECISE_SEAM_RIGHT, ModelVolumeType::PRECISE_SEAM_BLOCKED);
    SeamFixture fixture(3);
    PrintObject *object = fixture.print.get_object(0);
    const auto &slicing = object->slicing_parameters();
    REQUIRE(slicing.raft_layers() > 0);
    REQUIRE(fixture.layer->id() == slicing.raft_layers());
    // The real modifier mesh intersects the lower sampled object layer, but not the upper one.
    Layer *upper = object->add_layer(int(slicing.raft_layers() + 1), 0.2, slicing.object_print_z_min + 4.2, 4.1);
    ModelVolume *modifier = fixture.modifiers.objects.front()->add_volume(Test::cube(4));
    modifier->set_type(mode);
    const auto &slices = fixture.cache.emplace(modifier, PreciseSeam::prepare_modifier_slices(object->slice_single_volume_regions(modifier))).first->second;
    REQUIRE(slices.size() == 2);
    REQUIRE(slices[0].size() == 1);
    CHECK(slices[1].empty());

    // Position the test perimeter relative to the transformed slice to isolate layer indexing.
    const BoundingBox &bounds = slices[0].front().bounds;
    const Point origin = bounds.min;
    const Polygon original(Points{origin + mm(-2, 2), origin + mm(6, 2),
                                  origin + mm(6, 10), origin + mm(-2, 10)});
    Polygon perimeter = original;
    if (is_precise_seam_strong(mode)) {
        const auto seam = PreciseSeam::insert_strong_seam_point(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
        REQUIRE(seam.has_value());
        const double x = mode == ModelVolumeType::PRECISE_SEAM_LEFT ? 0. :
                         mode == ModelVolumeType::PRECISE_SEAM_RIGHT ? 4. : 2.;
        CHECK(*seam == origin + mm(x, 2));
        perimeter = original;
        CHECK_FALSE(PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), upper, fixture.cache).has_value());
    } else {
        const auto segments = PreciseSeam::collect_weak_modifier_segments(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
        REQUIRE(segments.size() == 1);
        CHECK(segments.front().left_point == origin + mm(0, 2));
        CHECK(segments.front().right_point == origin + mm(4, 2));
        const auto types = weak_candidate_types(perimeter, segments);
        for (size_t i = 0; i < perimeter.size(); ++i) {
            const Point local = perimeter[i] - origin;
            const bool inside = local.y() == scale_(2.) && local.x() >= 0 && local.x() <= scale_(4.);
            CHECK(types[i] == (inside ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                                       SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral));
        }
        perimeter = original;
        CHECK(PreciseSeam::collect_weak_modifier_segments({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), upper, fixture.cache).empty());
    }
    CHECK(perimeter.points == original.points);
}

TEST_CASE("Strong seam modes select the requested location on a clipped side", "[PreciseSeam]")
{
    const auto mode = GENERATE(ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_CENTER,
                              ModelVolumeType::PRECISE_SEAM_RIGHT);
    SeamFixture fixture;
    // Center must retain the correct source edge across several collinear edges.
    Polygon perimeter(Points{mm(0, 0), mm(2, 0), mm(4, 0), mm(8, 0), mm(20, 0), mm(20, 20), mm(0, 20)});
    const ModelVolume *modifier = fixture.add(mode, {rectangle(1, -2, 13, 2)});
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(seam.has_value());
    const double expected_x = mode == ModelVolumeType::PRECISE_SEAM_LEFT ? 1.0 :
                              mode == ModelVolumeType::PRECISE_SEAM_RIGHT ? 13.0 : 7.0;
    CHECK(*seam == mm(expected_x, 0));
    require_vertex(perimeter, mm(expected_x - 0.001, 0));
    require_vertex(perimeter, mm(expected_x + 0.001, 0));
    check_square_boundary(perimeter);
    CHECK_FALSE(warnings.full_containment.load());
    CHECK_FALSE(warnings.multiple_intersections.load());
}

TEST_CASE("Center seams preserve the edge order at vertices and across the contour origin", "[PreciseSeam]")
{
    const bool wrap = GENERATE(false, true);
    SeamFixture fixture;
    Polygon perimeter(Points{mm(0, 0), mm(4, 0), mm(20, 0), mm(20, 20), mm(0, 20)});
    // Symmetric cuts put the arc midpoint exactly on an existing vertex, including vertex zero.
    const Polygon cut = wrap ? rectangle(-2, -2, 4, 4) : rectangle(1, -2, 7, 2);
    const ModelVolume *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {cut});
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(seam.has_value());
    CHECK(*seam == (wrap ? mm(0, 0) : mm(4, 0)));
    require_vertex(perimeter, wrap ? mm(0, 0.001) : mm(3.999, 0));
    require_vertex(perimeter, wrap ? mm(0.001, 0) : mm(4.001, 0));
    check_square_boundary(perimeter);
}

TEST_CASE("Center seams land on the closing edge", "[PreciseSeam]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const auto *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(-2, 3, 2, 9)});
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(seam.has_value());
    CHECK(*seam == mm(0, 6));
    require_vertex(perimeter, mm(0, 5.999));
    require_vertex(perimeter, mm(0, 6.001));
    check_square_boundary(perimeter);
}

TEST_CASE("Weak boundaries on the closing edge preserve candidate types", "[PreciseSeam][Regression]")
{
    const bool reverse = GENERATE(false, true);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    // Keep vertex zero fixed so the selected side remains the closing edge in both directions.
    if (reverse)
        std::reverse(perimeter.points.begin() + 1, perimeter.points.end());
    const Polygon cut = reverse ? rectangle(3, -2, 9, 2) : rectangle(-2, 3, 2, 9);
    const auto *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {cut});
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 1);
    const auto types = weak_candidate_types(perimeter, segments);
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &p = perimeter[i];
        const coord_t along = reverse ? p.x() : p.y();
        const coord_t across = reverse ? p.y() : p.x();
        const bool inside = across == 0 && along >= scale_(3.) && along <= scale_(9.);
        CHECK(types[i] == (inside ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                                   SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral));
    }
    check_square_boundary(perimeter);
}

TEST_CASE("Strong intersection warnings count joined segments across vertex zero", "[PreciseSeam][Regression]")
{
    const bool extra_segment = GENERATE(false, true);
    const auto mode = GENERATE(ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_CENTER,
                              ModelVolumeType::PRECISE_SEAM_RIGHT);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    ExPolygons regions{ExPolygon(rectangle(-2, -2, 4, 4))};
    if (extra_segment)
        regions.emplace_back(rectangle(8, -2, 12, 2));
    // The two fragments at vertex zero form one eight-millimeter segment, longer than the extra one.
    const auto extracted = PreciseSeam::extract_perimeter_segments(
        PreciseSeam::PreparedPerimeter(perimeter), PreciseSeam::prepare_modifier_regions(regions), mode);
    REQUIRE(extracted.segments.size() == (extra_segment ? 2 : 1));
    const auto *modifier = fixture.add_regions(mode, std::move(regions));
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto seam = PreciseSeam::insert_strong_seam_point(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(seam.has_value());
    CHECK(*seam == (mode == ModelVolumeType::PRECISE_SEAM_LEFT ? mm(0, 4) :
                   mode == ModelVolumeType::PRECISE_SEAM_RIGHT ? mm(4, 0) : mm(0, 0)));
    CHECK(warnings.multiple_intersections.load() == (extra_segment ? PreciseSeam::PreciseSeamWarnings::type_bit(mode) : 0u));
    CHECK(warnings.failed_fragments.load() == 0);
    CHECK_FALSE(warnings.full_containment.load());
    check_square_boundary(perimeter);
}

TEST_CASE("A collinear contour origin preserves strong targets and weak candidate types", "[PreciseSeam][Regression]")
{
    const bool reverse = GENERATE(false, true);
    const auto mode = GENERATE(ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_CENTER,
                              ModelVolumeType::PRECISE_SEAM_RIGHT, ModelVolumeType::PRECISE_SEAM_BLOCKED);
    SeamFixture fixture;
    Polygon perimeter(Points{mm(10, 0), mm(20, 0), mm(20, 20), mm(0, 20), mm(0, 0)});
    // Both clipping fragments have two points, with the artificial cut inside a straight side.
    if (reverse)
        std::reverse(perimeter.points.begin() + 1, perimeter.points.end());
    const auto *modifier = fixture.add(mode, {rectangle(8, -2, 12, 2)});
    PreciseSeam::PreciseSeamWarnings warnings;
    if (is_precise_seam_strong(mode)) {
        const auto seam = PreciseSeam::insert_strong_seam_point(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        REQUIRE(seam.has_value());
        const double x = mode == ModelVolumeType::PRECISE_SEAM_CENTER ? 10. :
            ((mode == ModelVolumeType::PRECISE_SEAM_LEFT) != reverse ? 8. : 12.);
        CHECK(*seam == mm(x, 0));
    } else {
        const auto segments = PreciseSeam::collect_weak_modifier_segments(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        REQUIRE(segments.size() == 1);
        const auto types = weak_candidate_types(perimeter, segments);
        for (size_t i = 0; i < perimeter.size(); ++i) {
            const Point &p = perimeter[i];
            const bool inside = p.y() == 0 && p.x() >= scale_(8.) && p.x() <= scale_(12.);
            CHECK(types[i] == (inside ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                                       SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral));
        }
    }
    CHECK_FALSE(warnings.multiple_intersections.load());
    CHECK_FALSE(warnings.full_containment.load());
    CHECK(warnings.failed_fragments.load() == 0);
    check_square_boundary(perimeter);
}

TEST_CASE("A gap on the closing edge preserves the complementary seam segment", "[PreciseSeam][Regression]")
{
    const auto mode = GENERATE(ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_CENTER,
                              ModelVolumeType::PRECISE_SEAM_RIGHT, ModelVolumeType::PRECISE_SEAM_BLOCKED);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    ExPolygon region(rectangle(-2, -2, 22, 22));
    region.holes.push_back(rectangle(-1, 8, 1, 12));
    region.holes.back().reverse();
    // The retained 76 mm arc starts and ends on edge n-1, with a four-millimeter gap between them.
    const auto *modifier = fixture.add_regions(mode, {region});
    PreciseSeam::PreciseSeamWarnings warnings;
    if (is_precise_seam_strong(mode)) {
        const auto seam = PreciseSeam::insert_strong_seam_point(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        REQUIRE(seam.has_value());
        CHECK(*seam == (mode == ModelVolumeType::PRECISE_SEAM_LEFT ? mm(0, 8) :
                       mode == ModelVolumeType::PRECISE_SEAM_RIGHT ? mm(0, 12) : mm(20, 10)));
    } else {
        const auto segments = PreciseSeam::collect_weak_modifier_segments(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        REQUIRE(segments.size() == 1);
        const auto types = weak_candidate_types(perimeter, segments);
        for (size_t i = 0; i < perimeter.size(); ++i) {
            const Point &p = perimeter[i];
            const bool gap = p.x() == 0 && p.y() > scale_(8.) && p.y() < scale_(12.);
            CHECK(types[i] == (gap ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral :
                                    SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked));
        }
    }
    CHECK_FALSE(warnings.multiple_intersections.load());
    CHECK_FALSE(warnings.full_containment.load());
    CHECK(warnings.failed_fragments.load() == 0);
    check_square_boundary(perimeter);
}

TEST_CASE("Enforced oversampling stays inside a zone that wraps around a gap", "[PreciseSeam][Regression]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    ExPolygon region(rectangle(-2, -2, 22, 22));
    region.holes.push_back(rectangle(-1, 8, 1, 12));
    region.holes.back().reverse();
    // The zone runs from (0, 8) around the square to (0, 12); the edge after its right boundary is the gap.
    // Boundary helpers already keep that edge too short to split, so this guards the invariant, not a visible bug.
    const auto *modifier = fixture.add_regions(ModelVolumeType::PRECISE_SEAM_ENFORCED, {region});
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 1);
    CHECK(segments[0].left_point == mm(0, 8));
    CHECK(segments[0].right_point == mm(0, 12));

    const auto types = weak_candidate_types(perimeter, segments);
    Points gap_points;
    size_t enforced_count = 0;
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &p = perimeter[i];
        const bool gap = p.x() == 0 && p.y() > scale_(8.) && p.y() < scale_(12.);
        if (gap)
            gap_points.push_back(p);
        else
            ++enforced_count;
        CHECK(types[i] == (gap ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral :
                                SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced));
    }
    // Only the two boundary refinement helpers may lie in the gap; no oversampling points.
    std::sort(gap_points.begin(), gap_points.end(), [](const Point &a, const Point &b) { return a.y() < b.y(); });
    CHECK(gap_points == Points{mm(0, 8.001), mm(0, 11.999)});
    // The 76 mm zone itself is oversampled: far more candidates than its corners and boundaries.
    CHECK(enforced_count > 76. / SeamPlacer::enforcer_oversampling_distance - 10);
    check_square_boundary(perimeter);
}

TEST_CASE("Coincident weak boundaries do not prevent later boundary refinement", "[PreciseSeam][Regression]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    // Duplicate boundaries must not stop refinement before the separate segment is reached.
    const auto *a = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(2, -2, 6, 2)});
    const auto *b = fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, {rectangle(2, -2, 6, 2)});
    const auto *c = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(10, -2, 14, 2)});
    const auto segments = PreciseSeam::collect_weak_modifier_segments({c, b, a}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 3);
    for (double x : {1.999, 6.001, 9.999, 14.001})
        require_vertex(perimeter, mm(x, 0));
    check_square_boundary(perimeter);
}

TEST_CASE("Weak boundaries on one source edge retain insertion order and modifier priority", "[PreciseSeam][Regression]")
{
    const bool closing_edge = GENERATE(false, true);
    const bool coincident = GENERATE(false, true);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Polygon outer = closing_edge ? rectangle(-2, 2, 2, 8) : rectangle(2, -2, 8, 2);
    const Polygon inner = coincident ? outer :
        (closing_edge ? rectangle(-2, 4, 2, 6) : rectangle(4, -2, 6, 2));
    const auto *blocked = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {outer});
    const auto *neutral = fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, {inner});
    // Priority order is not geometric insertion order; Neutral must win in the overlap.
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {blocked, neutral}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 2);
    CHECK(segments[0].type == SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked);
    CHECK(segments[1].type == SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral);
    CHECK(segments[0].left_position.edge_index == (closing_edge ? 3 : 0));
    CHECK(segments[0].right_position.edge_index == (closing_edge ? 3 : 0));
    CHECK_THAT(segments[0].left_position.parameter, Catch::Matchers::WithinAbs(closing_edge ? 0.6 : 0.1, 1e-12));
    CHECK_THAT(segments[0].right_position.parameter, Catch::Matchers::WithinAbs(closing_edge ? 0.9 : 0.4, 1e-12));

    const auto types = weak_candidate_types(perimeter, segments);
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &p = perimeter[i];
        const coord_t along = closing_edge ? p.y() : p.x();
        const coord_t across = closing_edge ? p.x() : p.y();
        const bool in_outer = across == 0 && along >= scale_(2.) && along <= scale_(8.);
        const bool in_inner = coincident ? in_outer :
            (across == 0 && along >= scale_(4.) && along <= scale_(6.));
        const auto expected = in_outer && !in_inner ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                                                     SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral;
        CHECK(types[i] == expected);
    }
    // Ascending insertion on one edge would retrace the contour or lose pending boundaries.
    check_square_boundary(perimeter);
}

TEST_CASE("Weak zones overwrite painted candidate types inside their boundaries only", "[PreciseSeam][Regression]")
{
    using Type = SeamPlacerImpl::EnforcedBlockedSeamPoint;
    // The bottom side is painted Enforced or Blocked, the top side Blocked; the vertical sides stay Neutral.
    const Type painted_bottom = GENERATE(Type::Enforced, Type::Blocked);
    const auto zone_type = GENERATE(ModelVolumeType::PRECISE_SEAM_ENFORCED, ModelVolumeType::PRECISE_SEAM_BLOCKED,
                                    ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    CAPTURE(int(painted_bottom), int(zone_type));
    const Type expected_zone = zone_type == ModelVolumeType::PRECISE_SEAM_ENFORCED ? Type::Enforced :
                               zone_type == ModelVolumeType::PRECISE_SEAM_BLOCKED  ? Type::Blocked : Type::Neutral;
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    // The zone covers x in [6, 14] on the bottom side only.
    const auto *modifier = fixture.add(zone_type, {rectangle(6, -2, 14, 2)});
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 1);
    const auto painted = [painted_bottom](const Point &p) {
        return p.y() == 0 ? painted_bottom : p.y() == mm(0, 20).y() ? Type::Blocked : Type::Neutral;
    };
    const auto types = weak_candidate_types(perimeter, segments, painted);
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &p = perimeter[i];
        CAPTURE(p.x(), p.y());
        const bool inside = p.y() == 0 && p.x() >= mm(6, 0).x() && p.x() <= mm(14, 0).x();
        // Inside the zone the weak type wins over painting, Neutral clearing it; outside, painting stays.
        CHECK(types[i] == (inside ? expected_zone : painted(p)));
    }
}

TEST_CASE("Weak boundaries sharing a vertex refine both sides", "[PreciseSeam]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const auto *a = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(2, -2, 6, 2)});
    const auto *b = fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, {rectangle(6, -2, 10, 2)});
    const auto segments = PreciseSeam::collect_weak_modifier_segments({a, b}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 2);
    require_vertex(perimeter, mm(5.999, 0));
    require_vertex(perimeter, mm(6.001, 0));
    check_square_boundary(perimeter);
}

TEST_CASE("Weak processing applies every ready interval and leaves gaps unchanged", "[PreciseSeam][SegmentExtraction]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    ExPolygon area(rectangle(1, -3, 10, 3));
    area.holes.push_back(rectangle(4, -1, 7, 1));
    area.holes.back().reverse();
    const auto *modifier = fixture.add_regions(ModelVolumeType::PRECISE_SEAM_BLOCKED,
                                               {area, ExPolygon(rectangle(14, -3, 18, 3))});
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(segments.size() == 3);
    const Points starts{mm(1, 0), mm(7, 0), mm(14, 0)};
    const Points ends{mm(4, 0), mm(10, 0), mm(18, 0)};
    for (size_t i = 0; i < segments.size(); ++i) {
        CHECK(segments[i].left_point == starts[i]);
        CHECK(segments[i].right_point == ends[i]);
        require_vertex(perimeter, starts[i]);
        require_vertex(perimeter, ends[i]);
    }
    const auto types = weak_candidate_types(perimeter, segments);
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &point = perimeter[i];
        // Classify the three independent intervals, including their boundary vertices.
        const bool inside = point.y() == 0 &&
            ((point.x() >= starts[0].x() && point.x() <= ends[0].x()) ||
             (point.x() >= starts[1].x() && point.x() <= ends[1].x()) ||
             (point.x() >= starts[2].x() && point.x() <= ends[2].x()));
        CHECK(types[i] == (inside ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                                  SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral));
    }
    CHECK_FALSE(warnings.multiple_intersections.load());
    CHECK_FALSE(warnings.full_containment.load());
    CHECK(warnings.failed_fragments.load() == 0);
    check_square_boundary(perimeter);
}

TEST_CASE("Weak priority and enforcement refinement apply on both crossed sides", "[PreciseSeam][SegmentExtraction]")
{
    const auto high_type = GENERATE(ModelVolumeType::PRECISE_SEAM_BLOCKED, ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    const auto expected_high = high_type == ModelVolumeType::PRECISE_SEAM_BLOCKED ?
        SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked : SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral;
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const auto *low = fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, {rectangle(8, -2, 12, 22)});
    const auto *high = fixture.add(high_type, {rectangle(10, -2, 14, 22)});
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {low, high}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 4);
    CHECK(segments[0].type == SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced);
    CHECK(segments[1].type == SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced);
    CHECK(segments[2].type == expected_high);
    CHECK(segments[3].type == expected_high);
    const auto types = weak_candidate_types(perimeter, segments);
    size_t enforced_counts[2] = {0, 0};
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &point = perimeter[i];
        const bool side = point.y() == 0 || point.y() == mm(0, 20).y();
        auto expected = SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral;
        if (side && point.x() >= mm(8, 0).x() && point.x() <= mm(12, 0).x())
            expected = SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced;
        if (side && point.x() >= mm(10, 0).x() && point.x() <= mm(14, 0).x())
            expected = expected_high;
        CHECK(types[i] == expected);
        if (expected == SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced)
            ++enforced_counts[point.y() == 0 ? 0 : 1];
    }
    // Both surviving enforced patches must contain interior candidates, not just boundaries.
    CHECK(enforced_counts[0] > 2);
    CHECK(enforced_counts[1] > 2);
    check_square_boundary(perimeter);
}

TEST_CASE("Weak full containment follows painting: the whole perimeter takes the zone's type", "[PreciseSeam]")
{
    const auto type = GENERATE(ModelVolumeType::PRECISE_SEAM_ENFORCED, ModelVolumeType::PRECISE_SEAM_BLOCKED,
                               ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    CAPTURE(type);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Points original = perimeter.points;
    const auto *modifier = fixture.add(type, {rectangle(-2, -2, 22, 22)});
    PreciseSeam::PreciseSeamWarnings warnings;
    // Registered as SeamPlacer::init() does, so the per-modifier warning flag is tracked.
    warnings.modifier_usage.try_emplace(modifier);
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    // Painting of both kinds, so that each type's effect on it is visible.
    const auto painted = [](const Point &p) {
        return p.x() < scale_(10.) ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced :
                                     SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked;
    };
    const auto types = weak_candidate_types(perimeter, segments, painted, &warnings);
    CHECK_FALSE(warnings.multiple_intersections.load());
    CHECK(warnings.failed_fragments.load() == 0);
    REQUIRE(segments.size() == 1);
    CHECK(segments.front().whole_perimeter);
    CHECK(warnings.full_containment.load() == 0);
    // Only a whole Blocked perimeter is reported, as the seam cannot avoid it.
    CHECK(warnings.modifier_usage.at(modifier).fully_blocked.load() == (type == ModelVolumeType::PRECISE_SEAM_BLOCKED));
    // Every candidate takes the zone's type, overriding painting.
    const auto expected = type == ModelVolumeType::PRECISE_SEAM_ENFORCED ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced :
                          type == ModelVolumeType::PRECISE_SEAM_BLOCKED  ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                                                                          SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral;
    for (size_t i = 0; i < perimeter.size(); ++i)
        CHECK(types[i] == expected);
    if (type == ModelVolumeType::PRECISE_SEAM_ENFORCED) {
        // As if painted green all round: every 20 mm edge is subdivided into steps of at most 0.2 mm.
        CHECK(perimeter.size() >= size_t(80.f / SeamPlacer::enforcer_oversampling_distance));
        check_square_boundary(perimeter);
    } else
        CHECK(perimeter.points == original); // Like a perimeter unmarked or painted red all round: no subdivision.
}

TEST_CASE("Whole-perimeter weak zones take part in the usual priority order", "[PreciseSeam]")
{
    // Volumes are listed low priority first; a later zone overwrites an earlier one.
    const int scenario = GENERATE(0, 1, 2, 3, 4, 5, 6);
    CAPTURE(scenario);
    SeamFixture fixture;
    const auto *whole_enforced = fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, {rectangle(-2, -2, 22, 22)});
    const auto *whole_neutral = fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, {rectangle(-2, -2, 22, 22)});
    const auto *whole_blocked = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(-2, -2, 22, 22)});
    const auto *lower_whole_blocked = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(-2, -2, 22, 22)});
    const auto *band_blocked = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(8, -2, 12, 2)});
    const auto *band_enforced = fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, {rectangle(8, -2, 12, 2)});
    std::vector<const ModelVolume*> volumes;
    if (scenario == 0) volumes = {whole_enforced, band_blocked}; // A band forbidden inside a green perimeter.
    if (scenario == 1) volumes = {band_enforced, whole_neutral}; // A whole Neutral clears the lower band.
    if (scenario == 2) volumes = {band_enforced, whole_blocked}; // A whole Blocked overrides the lower band.
    if (scenario == 3) volumes = {whole_blocked, band_enforced}; // A higher band unblocks part of a whole Blocked.
    if (scenario == 4) volumes = {whole_blocked, band_blocked};  // A higher Blocked band keeps it whole Blocked.
    if (scenario == 5) volumes = {lower_whole_blocked, whole_enforced, whole_blocked}; // The lower one is overwritten.
    if (scenario == 6) volumes = {lower_whole_blocked, whole_blocked};                 // Only the higher one is in effect.
    Polygon perimeter = rectangle(0, 0, 20, 20);
    PreciseSeam::PreciseSeamWarnings warnings;
    // Registered as SeamPlacer::init() does, so the per-modifier warning flag is tracked.
    warnings.modifier_usage.try_emplace(whole_blocked);
    warnings.modifier_usage.try_emplace(lower_whole_blocked);
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        volumes, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    const auto types = weak_candidate_types(perimeter, segments, {}, &warnings);
    using Type = SeamPlacerImpl::EnforcedBlockedSeamPoint;
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Point &p = perimeter[i];
        CAPTURE(p.x(), p.y());
        const bool in_band = p.y() == 0 && p.x() >= scale_(8.) && p.x() <= scale_(12.);
        const Type expected = scenario == 0 ? (in_band ? Type::Blocked : Type::Enforced) :
                              scenario == 1 ? Type::Neutral :
                              scenario == 3 ? (in_band ? Type::Enforced : Type::Blocked) :
                                              Type::Blocked;
        CHECK(types[i] == expected);
    }
    CHECK(warnings.full_containment.load() == 0);
    // The warning follows the final typing: only a perimeter left whole Blocked is reported, through the
    // highest whole zone.
    const bool expected_fully_blocked = scenario == 2 || scenario == 4 || scenario == 5 || scenario == 6;
    CHECK(warnings.modifier_usage.at(whole_blocked).fully_blocked.load() == expected_fully_blocked);
    CHECK_FALSE(warnings.modifier_usage.at(lower_whole_blocked).fully_blocked.load());
    check_square_boundary(perimeter);
}

TEST_CASE("A partial Blocked zone on a perimeter painted red all round gives no whole-perimeter warning", "[PreciseSeam]")
{
    SeamFixture fixture;
    const auto *band_blocked = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(8, -2, 12, 2)});
    Polygon perimeter = rectangle(0, 0, 20, 20);
    PreciseSeam::PreciseSeamWarnings warnings;
    warnings.modifier_usage.try_emplace(band_blocked);
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {band_blocked}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(segments.size() == 1);
    // Red painting makes every candidate Blocked, so the check reaches the whole-zone filter.
    const auto painted = [](const Point &) { return SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked; };
    const auto types = weak_candidate_types(perimeter, segments, painted, &warnings);
    for (const auto type : types)
        CHECK(type == SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked);
    // Only a modifier that covers the whole perimeter is reported.
    CHECK_FALSE(warnings.modifier_usage.at(band_blocked).fully_blocked.load());
}

TEST_CASE("Coverage short by a gap under 1 um is full containment while a narrow band keeps its point", "[PreciseSeam][Regression]")
{
    // A tip pokes through the bottom edge at x = 10 mm by `depth` nm, so the perimeter line crosses it
    // over about 2 * depth. As a hole in a large region it leaves that much uncovered (a touch from
    // inside); as a separate triangle it is that much covered (a narrow band or an outside contact).
    const auto type = GENERATE(ModelVolumeType::PRECISE_SEAM_ENFORCED, ModelVolumeType::PRECISE_SEAM_CENTER,
                               ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_RIGHT);
    const bool gap = GENERATE(true, false);
    const coord_t depth = GENERATE(coord_t(2), coord_t(2000));
    CAPTURE(type, gap, depth);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Points original = perimeter.points;
    const Point tip(mm(10, 0).x(), gap ? -depth : depth);
    ExPolygon region;
    if (gap) {
        region = ExPolygon(rectangle(-2, -2, 22, 22));
        region.holes.push_back(Polygon(Points{tip, mm(11, 1), mm(9, 1)}));
        region.holes.back().reverse();
    } else {
        region = ExPolygon(Polygon(Points{mm(9, -1), mm(11, -1), tip}));
    }
    const auto *modifier = fixture.add_regions(type, {region});
    const bool micro = depth < scale_(0.001);
    PreciseSeam::PreciseSeamWarnings warnings;
    if (is_precise_seam_strong(type)) {
        const auto seam = PreciseSeam::insert_strong_seam_point(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        if (gap && micro) {
            // Skipped like an exact touch, instead of a seam at the touch or on the opposite side.
            CHECK_FALSE(seam.has_value());
            CHECK(perimeter.points == original);
        } else {
            REQUIRE(seam.has_value());
            if (!gap) // The band is processed normally: the seam lands on it.
                CHECK((*seam - mm(10, 0)).cast<double>().norm() <= double(depth) + 1.);
        }
    } else {
        const auto segments = PreciseSeam::collect_weak_modifier_segments(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        const auto types = weak_candidate_types(perimeter, segments);
        const size_t enforced = size_t(std::count(types.begin(), types.end(), SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced));
        if (gap && micro) {
            // Full containment like an exact touch, instead of collapsing the intended zone into one forced
            // point: Enforced then types the whole perimeter, subdivided as if painted green all round.
            REQUIRE(segments.size() == 1);
            CHECK(segments.front().whole_perimeter);
            CHECK(enforced == types.size());
            CHECK(enforced > 300);
        } else if (gap) {
            REQUIRE(segments.size() == 1);
            CHECK(enforced > 300); // A 4 um gap is real: nearly the whole 80 mm perimeter, oversampled.
        } else {
            // The band is processed normally; below 1 um its boundaries snap into a single candidate.
            REQUIRE(segments.size() == 1);
            CHECK(enforced == (micro ? 1u : 2u));
        }
    }
    // Only the sub-micron gap is full containment, and only strong types skip it with a warning;
    // nothing here is a binding failure.
    CHECK((warnings.full_containment.load() != 0) == (gap && micro && is_precise_seam_strong(type)));
    CHECK(warnings.failed_fragments.load() == 0);
    check_square_boundary(perimeter);
}

TEST_CASE("A modifier touching the perimeter at one point keeps the full containment policy", "[PreciseSeam][Regression]")
{
    const auto type = GENERATE(ModelVolumeType::PRECISE_SEAM_CENTER, ModelVolumeType::PRECISE_SEAM_LEFT,
                               ModelVolumeType::PRECISE_SEAM_RIGHT, ModelVolumeType::PRECISE_SEAM_ENFORCED,
                               ModelVolumeType::PRECISE_SEAM_BLOCKED, ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Points original = perimeter.points;
    // The hole touches the perimeter only at (10, 0): no gap, so the contact must not create segments.
    ExPolygon region(rectangle(-2, -2, 22, 22));
    region.holes.push_back(Polygon(Points{mm(10, 0), mm(12, 2), mm(10, 4), mm(8, 2)}));
    region.holes.back().reverse();
    const auto *modifier = fixture.add_regions(type, {region});
    PreciseSeam::PreciseSeamWarnings warnings;
    if (is_precise_seam_strong(type))
        CHECK_FALSE(PreciseSeam::insert_strong_seam_point(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings).has_value());
    else {
        const auto segments = PreciseSeam::collect_weak_modifier_segments(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        // A weak modifier gives its type to the whole perimeter.
        REQUIRE(segments.size() == 1);
        CHECK(segments.front().whole_perimeter);
    }
    // Only a whole-perimeter Enforced zone subdivides the edges.
    if (type != ModelVolumeType::PRECISE_SEAM_ENFORCED)
        CHECK(perimeter.points == original);
    // Only strong types are skipped with the full-containment warning.
    CHECK(warnings.full_containment.load() == (is_precise_seam_strong(type) ? PreciseSeam::PreciseSeamWarnings::type_bit(type) : 0u));
    CHECK_FALSE(warnings.multiple_intersections.load());
    CHECK(warnings.failed_fragments.load() == 0);
}

TEST_CASE("Unsupported strong modifier sections are skipped with the appropriate warning", "[PreciseSeam]")
{
    const int scenario = GENERATE(0, 1, 2, 3);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Points original = perimeter.points;
    Polygons slices;
    if (scenario == 0) slices = {rectangle(30, 30, 40, 40)}; // Disjoint bounds.
    if (scenario == 1) slices = {rectangle(2, 2, 4, 4)};     // Wholly inside; no common boundary.
    if (scenario == 2) slices = {rectangle(-2, -2, 22, 22)}; // Contains the entire perimeter.
    if (scenario == 3) {
        Polygon hole = rectangle(2, 2, 4, 4);
        hole.reverse();
        slices = {rectangle(-2, -2, 22, 22), hole};
    }
    const auto *modifier = scenario == 3 ?
        fixture.add_regions(ModelVolumeType::PRECISE_SEAM_CENTER, {ExPolygon(slices[0], slices[1])}) :
        fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, std::move(slices));
    PreciseSeam::PreciseSeamWarnings warnings;
    CHECK_FALSE(PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings).has_value());
    CHECK(perimeter.points == original);
    CHECK((warnings.full_containment.load() != 0) == (scenario == 2 || scenario == 3));
    CHECK_FALSE(warnings.multiple_intersections.load());
}

TEST_CASE("Strong selects the longest arc rather than the longest chord", "[PreciseSeam][SegmentExtraction]")
{
    const auto mode = GENERATE(ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_CENTER,
                               ModelVolumeType::PRECISE_SEAM_RIGHT);
    const size_t origin = GENERATE(size_t(0), size_t(2));
    const bool reverse_regions = GENERATE(false, true);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    std::rotate(perimeter.points.begin(), perimeter.points.begin() + origin, perimeter.points.end());
    // The corner arc is 8 mm with a shorter chord than the other 6 mm interval.
    ExPolygons regions{ExPolygon(rectangle(10, -2, 16, 2)), ExPolygon(rectangle(-2, -2, 4, 4))};
    if (reverse_regions) std::reverse(regions.begin(), regions.end());
    const auto *modifier = fixture.add_regions(mode, std::move(regions));
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(seam.has_value());
    CHECK(*seam == (mode == ModelVolumeType::PRECISE_SEAM_LEFT ? mm(0, 4) :
                   mode == ModelVolumeType::PRECISE_SEAM_RIGHT ? mm(4, 0) : mm(0, 0)));
    CHECK(warnings.multiple_intersections.load());
    check_square_boundary(perimeter);
}

TEST_CASE("Equal strong lengths compare the selected mode points", "[PreciseSeam][SegmentExtraction]")
{
    const auto mode = GENERATE(ModelVolumeType::PRECISE_SEAM_LEFT, ModelVolumeType::PRECISE_SEAM_CENTER,
                               ModelVolumeType::PRECISE_SEAM_RIGHT);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    // Both arcs are 6 mm: Left/Center prefer the corner; Right prefers the rear straight segment.
    const auto *modifier = fixture.add(mode, {rectangle(10, 18, 16, 22), rectangle(-2, 17, 3, 22)});
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(seam.has_value());
    CHECK(*seam == (mode == ModelVolumeType::PRECISE_SEAM_LEFT ? mm(3, 20) :
                   mode == ModelVolumeType::PRECISE_SEAM_RIGHT ? mm(10, 20) : mm(0, 20)));
    check_square_boundary(perimeter);
}

TEST_CASE("Strong length comparison does not merge nearly equal lengths", "[PreciseSeam][SegmentExtraction]")
{
    const coord_t extra = GENERATE(coord_t(0), coord_t(1));
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    Polygon second = rectangle(12, -2, 16, 2);
    second.points[1].x() += extra;
    second.points[2].x() += extra;
    const auto *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_LEFT, {rectangle(2, -2, 6, 2), second});
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(seam.has_value());
    CHECK(*seam == (extra == 0 ? mm(2, 0) : mm(12, 0)));
}

TEST_CASE("Strong priority wins over a longer segment in another modifier", "[PreciseSeam][SegmentExtraction]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const auto *high = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(2, -2, 4, 2)});
    const auto *low = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(1, 18, 19, 22)});
    const auto seam = PreciseSeam::insert_strong_seam_point({high, low}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(seam.has_value());
    CHECK(*seam == mm(3, 0));
    CHECK(std::find(perimeter.points.begin(), perimeter.points.end(), mm(10, 20)) == perimeter.points.end());
}

TEST_CASE("Strong continues past modifiers without usable segments", "[PreciseSeam][SegmentExtraction]")
{
    const bool contained = GENERATE(false, true);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const auto *first = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER,
        {contained ? rectangle(-2, -2, 22, 22) : rectangle(30, 30, 40, 40)});
    const auto *second = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(2, -2, 6, 2)});
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto seam = PreciseSeam::insert_strong_seam_point({first, second}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(seam.has_value());
    CHECK(*seam == mm(4, 0));
    CHECK((warnings.full_containment.load() != 0) == contained);
}

TEST_CASE("Strong compares all segments of structured modifier regions", "[PreciseSeam][SegmentExtraction]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    ExPolygon area(rectangle(1, -3, 10, 3));
    area.holes.push_back(rectangle(4, -1, 7, 1));
    area.holes.back().reverse();
    const auto *modifier = fixture.add_regions(ModelVolumeType::PRECISE_SEAM_CENTER,
                                               {area, ExPolygon(rectangle(14, -3, 18, 3))});
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(seam.has_value());
    CHECK(*seam == mm(16, 0));
    CHECK(warnings.multiple_intersections.load());
    CHECK_FALSE(warnings.full_containment.load());
    CHECK(warnings.failed_fragments.load() == 0);
    check_square_boundary(perimeter);
}

TEST_CASE("Strong tie order uses bed axes after rotating and translating the input", "[PreciseSeam][SegmentExtraction]")
{
    const bool rotated = GENERATE(false, true);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    Polygon modifier_region = rectangle(8, -2, 12, 22);
    const auto transform = [rotated](Point &p) {
        // Rotation is already baked into slice coordinates; translation cannot change ordering.
        if (rotated) p = Point(-p.y(), p.x());
        p += mm(100, 200);
    };
    for (Point &p : perimeter.points) transform(p);
    for (Point &p : modifier_region.points) transform(p);
    const auto *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {modifier_region});
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(seam.has_value());
    CHECK(*seam == (rotated ? mm(80, 210) : mm(110, 220)));
}

TEST_CASE("Strong modifiers warn when another slice polygon also intersects the perimeter", "[PreciseSeam]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    // Equal lengths on the same side select the leftmost mode point, independent of slice order.
    const auto *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER,
                                       {rectangle(2, -2, 6, 2), rectangle(12, -2, 16, 2)});
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    REQUIRE(seam.has_value());
    CHECK(*seam == mm(4, 0));
    CHECK(warnings.multiple_intersections.load());
    CHECK_FALSE(warnings.full_containment.load());
    check_square_boundary(perimeter);
}

TEST_CASE("Crossing modifiers choose the rear strong segment and keep both weak segments", "[PreciseSeam]")
{
    const bool strong = GENERATE(false, true);
    CAPTURE(strong);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    // The strip exits on opposite sides, leaving two exterior pieces. None of the clipped
    // vertices matches a square corner, so this also exercises the general segment-extraction path.
    const auto type = strong ? ModelVolumeType::PRECISE_SEAM_CENTER : ModelVolumeType::PRECISE_SEAM_BLOCKED;
    const auto *modifier = fixture.add(type, {rectangle(8, -2, 12, 22)});
    PreciseSeam::PreciseSeamWarnings warnings;
    if (strong) {
        const auto seam = PreciseSeam::insert_strong_seam_point({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        REQUIRE(seam.has_value());
        CHECK(*seam == mm(10, 20));
    } else {
        const auto segments = PreciseSeam::collect_weak_modifier_segments({modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
        REQUIRE(segments.size() == 2);
        CHECK(segments[0].left_point == mm(8, 0));
        CHECK(segments[0].right_point == mm(12, 0));
        CHECK(segments[1].left_point == mm(12, 20));
        CHECK(segments[1].right_point == mm(8, 20));
    }
    CHECK((warnings.multiple_intersections.load() != 0) == strong); // Warn only strong, after collecting all segments.
    CHECK_FALSE(warnings.full_containment.load());
    check_square_boundary(perimeter);
}

TEST_CASE("Modifier hierarchy keeps strong order and applies the highest weak priority last", "[PreciseSeam]")
{
    const auto high_type = GENERATE(ModelVolumeType::PRECISE_SEAM_BLOCKED, ModelVolumeType::PRECISE_SEAM_NEUTRAL,
                                   ModelVolumeType::PRECISE_SEAM_ENFORCED);
    const auto expected_type = high_type == ModelVolumeType::PRECISE_SEAM_BLOCKED ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Blocked :
                               high_type == ModelVolumeType::PRECISE_SEAM_NEUTRAL ? SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral :
                                                                                  SeamPlacerImpl::EnforcedBlockedSeamPoint::Enforced;
    SeamFixture fixture;
    const auto *strong_a = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(1, -2, 3, 2)});
    const auto *strong_b = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(11, -2, 13, 2)});
    const auto *high = fixture.add(high_type, {rectangle(2, -2, 6, 2)});
    const auto low_type = high_type == ModelVolumeType::PRECISE_SEAM_ENFORCED ? ModelVolumeType::PRECISE_SEAM_BLOCKED :
                                                                                           ModelVolumeType::PRECISE_SEAM_ENFORCED;
    const auto *low = fixture.add(low_type, {rectangle(2, -2, 6, 2)});
    std::vector<const ModelVolume*> strong, weak;
    bool has_strong = false;
    PreciseSeam::init_precise_seam_data(strong, weak, has_strong, fixture.modifiers.objects.front());
    REQUIRE(has_strong);
    CHECK(strong == std::vector<const ModelVolume*>{strong_a, strong_b});
    CHECK(weak == std::vector<const ModelVolume*>{low, high});
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const auto seam = PreciseSeam::insert_strong_seam_point(strong, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(seam.has_value());
    CHECK(*seam == mm(2, 0));

    perimeter = rectangle(0, 0, 20, 20);
    const auto segments = PreciseSeam::collect_weak_modifier_segments(weak, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 2);
    PrintObjectSeamData::LayerSeams result;
    result.perimeters.emplace_back();
    auto &loop = result.perimeters.back();
    // Include a preceding candidate to exercise nonzero global layer indices.
    result.points.emplace_back(Vec3f(-1, -1, 0), loop, 0, SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral);
    loop.start_index = 1;
    for (const Point &p : perimeter.points) {
        // Match production's double-to-float conversion: weak boundary lookup uses exact equality.
        const Vec2f position = unscale(p).cast<float>();
        result.points.emplace_back(Vec3f(position.x(), position.y(), 0), loop, 0,
                                   SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral);
    }
    loop.end_index = result.points.size();
    bool enforced = false;
    PreciseSeam::apply_weak_modifiers_to_perimeter(segments, result, loop, enforced);
    size_t patch_count = 0;
    for (size_t i = loop.start_index; i < loop.end_index; ++i) {
        const auto &candidate = result.points[i];
        // Axis-aligned input and interpolation keep y exactly zero; this classifies, rather than measures, the patch.
        const bool in_patch = candidate.position.y() == 0 && candidate.position.x() >= 2 && candidate.position.x() <= 6;
        CHECK(candidate.type == (in_patch ? expected_type :
                                          SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral));
        if (in_patch) ++patch_count;
    }
    CHECK(patch_count >= 2);
    if (high_type == ModelVolumeType::PRECISE_SEAM_ENFORCED) {
        // Four millimetres of enforcement must be subdivided, not just marked at its endpoints.
        CHECK(enforced);
        CHECK(patch_count >= size_t(4.0f / SeamPlacer::enforcer_oversampling_distance));
    }
    CHECK(result.points.front().type == SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral);
}

TEST_CASE("Volume sorting keeps strong before weak Precise Seam modifiers and the user order inside each group", "[PreciseSeam][Model]")
{
    Model model;
    ModelObject *object = model.add_object();
    // Enum order inside a group (Center < Right, Enforced < Neutral) is deliberately reversed here:
    // only the group may move a volume, the user's order inside it must survive.
    const auto add = [object](ModelVolumeType type) {
        ModelVolume *volume = object->add_volume(Test::cube(1));
        volume->set_type(type);
        return volume;
    };
    ModelVolume *part = add(ModelVolumeType::MODEL_PART);
    ModelVolume *modifier = add(ModelVolumeType::PARAMETER_MODIFIER);
    ModelVolume *second_part = add(ModelVolumeType::MODEL_PART);
    ModelVolume *weak_neutral = add(ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    ModelVolume *strong_right = add(ModelVolumeType::PRECISE_SEAM_RIGHT);
    ModelVolume *weak_enforced = add(ModelVolumeType::PRECISE_SEAM_ENFORCED);
    ModelVolume *strong_center = add(ModelVolumeType::PRECISE_SEAM_CENTER);
    const ModelVolumePtrs original = object->volumes;

    // A full sort orders parts before modifiers; a partial one keeps parts and modifiers as they are.
    const bool full_sort = GENERATE(false, true);
    CAPTURE(full_sort);
    object->volumes = original;
    object->sort_volumes(full_sort);
    const ModelVolumePtrs expected_head = full_sort ? ModelVolumePtrs{part, second_part, modifier}
                                                    : ModelVolumePtrs{part, modifier, second_part};
    ModelVolumePtrs expected = expected_head;
    for (ModelVolume *volume : {strong_right, strong_center, weak_neutral, weak_enforced})
        expected.push_back(volume);
    CHECK(object->volumes == expected);
}

TEST_CASE("Restoring Precise Seam positions overrides alignment only on perimeters with a strong point", "[PreciseSeam]")
{
    std::vector<PrintObjectSeamData::LayerSeams> layers(1);
    PrintObjectSeamData::LayerSeams &layer = layers.front();
    // Two perimeters of three candidates each; alignment has already finalized both.
    for (size_t loop_idx = 0; loop_idx < 2; ++loop_idx) {
        layer.perimeters.emplace_back();
        SeamPlacerImpl::Perimeter &loop = layer.perimeters.back();
        loop.start_index = layer.points.size();
        for (size_t i = 0; i < 3; ++i)
            layer.points.emplace_back(Vec3f(float(i), float(loop_idx), 0.f), loop, 0.f,
                                      SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral);
        loop.end_index = layer.points.size();
        loop.finalized = true;
        loop.seam_index = loop.start_index;
        loop.final_seam_position = Vec3f(0.5f, float(loop_idx), 0.f);
    }
    SeamPlacerImpl::Perimeter &strong = layer.perimeters[0];
    SeamPlacerImpl::Perimeter &plain = layer.perimeters[1];
    strong.precise_seam_point = Vec3f(2.f, 0.f, 0.f);
    strong.precise_seam_index = 2;

    PreciseSeam::restore_precise_seam_positions(layers);

    // The strong point and its candidate index replace the aligned position.
    CHECK(strong.final_seam_position == Vec3f(2.f, 0.f, 0.f));
    CHECK(strong.seam_index == 2);
    // A perimeter without a strong point keeps whatever alignment chose.
    CHECK(plain.final_seam_position == Vec3f(0.5f, 1.f, 0.f));
    CHECK(plain.seam_index == plain.start_index);
}

TEST_CASE("Modifier usage marks evaluated modifiers that never reach a perimeter", "[PreciseSeam]")
{
    SeamFixture fixture;
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto checked = [&warnings](const ModelVolume *volume) { return warnings.modifier_usage.at(volume).checked.load(); };
    const auto reached = [&warnings](const ModelVolume *volume) { return warnings.modifier_usage.at(volume).reached.load(); };

    // Weak: all modifiers are evaluated. One stops 1 mm short of the wall, one crosses it, one contains it.
    const ModelVolume *short_weak = fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, {rectangle(21, 5, 25, 10)});
    const ModelVolume *crossing_weak = fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, {rectangle(18, 5, 25, 10)});
    const ModelVolume *containing_weak = fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, {rectangle(-1, -1, 21, 21)});
    for (const ModelVolume *volume : {short_weak, crossing_weak, containing_weak})
        warnings.modifier_usage.try_emplace(volume);
    Polygon weak_perimeter = rectangle(0, 0, 20, 20);
    PreciseSeam::collect_weak_modifier_segments({short_weak, crossing_weak, containing_weak}, weak_perimeter,
        PreciseSeam::PreparedPerimeter(weak_perimeter), fixture.layer, fixture.cache, &warnings);
    CHECK(checked(short_weak));
    CHECK_FALSE(reached(short_weak));
    CHECK(reached(crossing_weak));
    // Full containment counts as reaching the perimeter.
    CHECK(reached(containing_weak));

    // Strong: evaluation stops at the first modifier with a segment. The one ending 0.5 mm before the
    // wall is evaluated without reaching it; the one after the winner is never evaluated, so nothing is
    // known about it and it must not be reported.
    const ModelVolume *short_strong = fixture.add(ModelVolumeType::PRECISE_SEAM_LEFT, {rectangle(5, -2, 10, -0.5)});
    const ModelVolume *winner = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(5, -2, 10, 2)});
    const ModelVolume *shadowed = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(12, -2, 15, 2)});
    for (const ModelVolume *volume : {short_strong, winner, shadowed})
        warnings.modifier_usage.try_emplace(volume);
    Polygon strong_perimeter = rectangle(0, 0, 20, 20);
    REQUIRE(PreciseSeam::insert_strong_seam_point({short_strong, winner, shadowed}, strong_perimeter,
        PreciseSeam::PreparedPerimeter(strong_perimeter), fixture.layer, fixture.cache, &warnings).has_value());
    CHECK(checked(short_strong));
    CHECK_FALSE(reached(short_strong));
    CHECK(reached(winner));
    CHECK_FALSE(checked(shadowed));

    // Unregistered modifiers are not tracked.
    const ModelVolume *unregistered = fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, {rectangle(18, 5, 25, 10)});
    Polygon other_perimeter = rectangle(0, 0, 20, 20);
    PreciseSeam::collect_weak_modifier_segments({unregistered}, other_perimeter,
        PreciseSeam::PreparedPerimeter(other_perimeter), fixture.layer, fixture.cache, &warnings);
    CHECK(warnings.modifier_usage.count(unregistered) == 0);
}

TEST_CASE("Weak zones type painting's oversampled candidates between their boundaries", "[PreciseSeam]")
{
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const auto *modifier = fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, {rectangle(8, -2, 12, 2)});
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache);
    REQUIRE(segments.size() == 1);

    // Candidates as production builds them for a bottom edge painted green: every polygon point, then
    // oversampled points every enforcer_oversampling_distance towards the next one, in float arithmetic.
    using Type = SeamPlacerImpl::EnforcedBlockedSeamPoint;
    PrintObjectSeamData::LayerSeams candidates;
    candidates.perimeters.emplace_back();
    auto &loop = candidates.perimeters.back();
    loop.start_index = 0;
    size_t oversampled = 0;
    for (size_t i = 0; i < perimeter.size(); ++i) {
        const Vec2f a = unscale(perimeter[i]).cast<float>();
        const Vec2f b = unscale(perimeter[(i + 1) % perimeter.size()]).cast<float>();
        const Vec3f position(a.x(), a.y(), 0.f);
        candidates.points.emplace_back(position, loop, 0.f, a.y() == 0.f ? Type::Enforced : Type::Neutral);
        if (a.y() != 0.f || b.y() != 0.f)
            continue;
        const Vec3f next(b.x(), b.y(), 0.f);
        const float distance = (next - position).norm();
        const Vec3f direction = (next - position).normalized();
        for (float step = SeamPlacer::enforcer_oversampling_distance; step < distance; step += SeamPlacer::enforcer_oversampling_distance) {
            candidates.points.emplace_back(position + direction * step, loop, 0.f, Type::Enforced);
            ++oversampled;
        }
    }
    loop.end_index = candidates.points.size();
    bool enforced = false;
    PreciseSeam::apply_weak_modifiers_to_perimeter(segments, candidates, loop, enforced);

    // Boundaries are found among the interleaved points by exact float identity; every candidate
    // between them, oversampled ones included, is cleared, and painting stays outside.
    size_t cleared_oversampled = 0;
    for (const auto &candidate : candidates.points) {
        const Vec3f &p = candidate.position;
        CAPTURE(p.x(), p.y());
        const bool in_zone = p.y() == 0.f && p.x() >= 8.f && p.x() <= 12.f;
        const Type expected = in_zone ? Type::Neutral : p.y() == 0.f ? Type::Enforced : Type::Neutral;
        CHECK(candidate.type == expected);
        if (in_zone && p.x() > 8.001f && p.x() < 11.999f)
            ++cleared_oversampled;
    }
    CHECK(oversampled > 0);
    CHECK(cleared_oversampled >= 15); // About 4 mm of 0.2 mm steps inside the zone.
}

namespace {
// A 20 x 20 x 2 mm cube sliced at 0.2 mm, so ten layers, each with one outer wall loop. Another object mesh
// takes the place of the cube.
struct ExportFixture {
    Model model;
    Print print;
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    Vec3d cube_min;

    explicit ExportFixture(const std::string &seam_position, int raft_layers = 0, TriangleMesh object = make_cube(20, 20, 2))
    {
        config.set_deserialize_strict("seam_position", seam_position);
        config.set_deserialize_strict("seam_slope_type", "none"); // The loop then starts exactly at the seam.
        config.set_deserialize_strict("layer_height", "0.2");
        config.set_deserialize_strict("initial_layer_print_height", "0.2");
        config.set_deserialize_strict("outer_wall_line_width", "0.4");
        config.set_deserialize_strict("initial_layer_line_width", "0.4");
        config.set_deserialize_strict("raft_layers", std::to_string(raft_layers));
        Test::init_print({std::move(object)}, print, model, config);
        const ModelVolume *cube = model.objects.front()->volumes.front();
        cube_min = cube->mesh().transformed_bounding_box(cube->get_matrix()).min;
    }

    // Places the helper's minimum corner at `min`, relative to the cube's minimum corner.
    void add(ModelVolumeType type, TriangleMesh mesh, const Vec3d &min, const std::string &name = "helper")
    {
        ModelVolume *volume = model.objects.front()->add_volume(std::move(mesh));
        volume->set_type(type);
        volume->name = name;
        const Vec3d current = volume->mesh().transformed_bounding_box(volume->get_matrix()).min;
        volume->set_offset(volume->get_offset() + cube_min + min - current);
    }

    std::string export_gcode()
    {
        print.apply(model, config);
        return Test::gcode(print);
    }

    std::string warning() const
    {
        std::string text;
        for (const auto &warning : print.step_state_with_warnings(psGCodeExport).warnings)
            text += warning.message;
        return text;
    }
};

struct OuterLoop {
    Vec2d seam; // Relative to the loop's centre, which is the cube's centre.
    BoundingBoxf bounds;
};

std::vector<OuterLoop> outer_wall_loops(const std::string &gcode)
{
    std::vector<OuterLoop> loops;
    bool outer_wall = false;
    bool in_loop = false;
    GCodeReader parser;
    parser.parse_buffer(gcode, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        const std::string_view comment = line.comment();
        if (comment.find("TYPE:") != std::string_view::npos) {
            outer_wall = comment.find("Outer wall") != std::string_view::npos;
            in_loop = false;
        }
        if (!outer_wall || !line.extruding(self) || line.dist_XY(self) == 0.f)
            return;
        if (!in_loop) {
            loops.push_back({Vec2d(self.x(), self.y()), {}});
            loops.back().bounds.merge(loops.back().seam);
            in_loop = true;
        }
        loops.back().bounds.merge(Vec2d(line.new_X(self), line.new_Y(self)));
    });
    for (OuterLoop &loop : loops)
        loop.seam -= loop.bounds.center();
    return loops;
}

TriangleMesh helper_with_hole()
{
    // 8 x 6 x 4 mm with a 4 x 3 mm hole through every layer of the cube.
    TriangleMesh helper = make_cube(8, 6, 4);
    TriangleMesh cavity = make_cube(4, 3, 3);
    cavity.translate(2., 1.5, 0.5);
    cavity.flip_triangles();
    helper.merge(cavity);
    return helper;
}
} // namespace

TEST_CASE("Strong modifiers place the exported seam on the side they cross", "[PreciseSeam]")
{
    // The helper crosses the front face (y = 0) at x = 4..6 mm. The loop runs counterclockwise, so on
    // the front face Left is the x = 4 mm edge and Right the x = 6 mm edge.
    const auto [type, x] = GENERATE(table<ModelVolumeType, double>({
        {ModelVolumeType::PRECISE_SEAM_CENTER, 5.},
        {ModelVolumeType::PRECISE_SEAM_LEFT, 4.},
        {ModelVolumeType::PRECISE_SEAM_RIGHT, 6.}}));
    const std::string seam_position = GENERATE(as<std::string>{}, "back", "aligned", "nearest");
    const int raft_layers = GENERATE(0, 2);
    CAPTURE(seam_position, raft_layers);
    ExportFixture fixture(seam_position, raft_layers);
    fixture.add(type, make_cube(2, 6, 4), Vec3d(4, -3, -1));
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10); // Raft layers have no outer wall.
    for (const OuterLoop &loop : loops) {
        CHECK_THAT(loop.seam.x(), Catch::Matchers::WithinAbs(x - 10., 0.01));
        CHECK_THAT(loop.seam.y(), Catch::Matchers::WithinAbs(loop.bounds.min.y() - loop.bounds.center().y(), 0.01));
    }
    CHECK(fixture.warning().empty());
}

TEST_CASE("Weak modifiers move the exported seam into or out of their zone", "[PreciseSeam]")
{
    // INVALID adds no helper: a back seam then lies on the rear face (y = 20 mm).
    const auto type = GENERATE(ModelVolumeType::INVALID, ModelVolumeType::PRECISE_SEAM_ENFORCED,
                               ModelVolumeType::PRECISE_SEAM_BLOCKED);
    ExportFixture fixture("back");
    if (type == ModelVolumeType::PRECISE_SEAM_ENFORCED)
        fixture.add(type, make_cube(2, 6, 4), Vec3d(4, -3, -1)); // Front face, x = 4..6 mm.
    if (type == ModelVolumeType::PRECISE_SEAM_BLOCKED)
        fixture.add(type, make_cube(30, 6, 4), Vec3d(-5, 17, -1)); // Everything from y = 17 mm back.
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (const OuterLoop &loop : loops) {
        CAPTURE(loop.seam.x(), loop.seam.y());
        if (type == ModelVolumeType::PRECISE_SEAM_ENFORCED) {
            CHECK_THAT(loop.seam.y(), Catch::Matchers::WithinAbs(loop.bounds.min.y() - loop.bounds.center().y(), 0.01));
            CHECK(loop.seam.x() >= -6.01);
            CHECK(loop.seam.x() <= -3.99);
        } else if (type == ModelVolumeType::PRECISE_SEAM_BLOCKED)
            CHECK(loop.seam.y() < 7.);
        else
            CHECK(loop.seam.y() > 7.);
    }
    CHECK(fixture.warning().empty());
}

TEST_CASE("A higher-priority Enforced band suppresses the whole-perimeter Blocked export warning", "[PreciseSeam][Regression]")
{
    ExportFixture fixture("back");
    // The earlier model volume has higher priority; weak application reverses this order.
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front enforcer");
    fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, make_cube(30, 30, 4), Vec3d(-5, -5, -1), "whole blocker");
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (const OuterLoop &loop : loops) {
        CAPTURE(loop.seam.x(), loop.seam.y());
        // The band wins on every layer, so the empty warning is not from modifiers that did nothing.
        CHECK_THAT(loop.seam.y(), Catch::Matchers::WithinAbs(loop.bounds.min.y() - loop.bounds.center().y(), 0.01));
        CHECK(loop.seam.x() >= -6.01);
        CHECK(loop.seam.x() <= -3.99);
    }
    CHECK(fixture.warning().empty());
}

TEST_CASE("A strong modifier with a hole uses the left of two equal crossings", "[PreciseSeam]")
{
    // The hole splits the front face crossing into x = 2..4 and x = 8..10 mm. The lengths tie, both
    // points lie equally far back, so the left one wins.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, helper_with_hole(), Vec3d(2, -3, -1));
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (const OuterLoop &loop : loops) {
        CHECK_THAT(loop.seam.x(), Catch::Matchers::WithinAbs(3. - 10., 0.01));
        CHECK_THAT(loop.seam.y(), Catch::Matchers::WithinAbs(loop.bounds.min.y() - loop.bounds.center().y(), 0.01));
    }
    // The multiple-intersections reason names the type.
    CHECK(fixture.warning().find("Seam Center") != std::string::npos);
}

TEST_CASE("The export warning names a modifier that misses the outer wall centreline", "[PreciseSeam]")
{
    // The helper reaches 0.1 mm into the cube, short of the centreline 0.2 mm inside a 0.4 mm outer wall.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, make_cube(2, 3.1, 4), Vec3d(4, -3, -1), "short helper");
    fixture.export_gcode();
    CHECK(fixture.warning().find("short helper") != std::string::npos);
}

TEST_CASE("Full containment warns for Blocked but not for Enforced", "[PreciseSeam]")
{
    const auto type = GENERATE(ModelVolumeType::PRECISE_SEAM_ENFORCED, ModelVolumeType::PRECISE_SEAM_BLOCKED);
    ExportFixture fixture("back");
    fixture.add(type, make_cube(30, 30, 4), Vec3d(-5, -5, -1), "whole helper");
    const auto loops = outer_wall_loops(fixture.export_gcode());
    // The whole perimeter takes the type either way, so the seam still follows the back preference.
    REQUIRE(loops.size() == 10);
    for (const OuterLoop &loop : loops)
        CHECK(loop.seam.y() > 7.);
    if (type == ModelVolumeType::PRECISE_SEAM_BLOCKED)
        CHECK(fixture.warning().find("Seam Blocked modifier \"whole helper\" of") != std::string::npos);
    else
        CHECK(fixture.warning().empty());
}

TEST_CASE("A perimeter blocked all round gets the seam of an unmarked one", "[PreciseSeam][Regression]")
{
    // Equal Blocked candidates compare like Neutral ones, so random and aligned seams still vary and align.
    const std::string seam_position = GENERATE(as<std::string>{}, "random", "aligned", "aligned_back", "back", "nearest");
    CAPTURE(seam_position);
    ExportFixture unmarked(seam_position);
    const auto expected = outer_wall_loops(unmarked.export_gcode());
    ExportFixture blocked(seam_position);
    blocked.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, make_cube(30, 30, 4), Vec3d(-5, -5, -1));
    const auto loops = outer_wall_loops(blocked.export_gcode());
    REQUIRE(loops.size() == expected.size());
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i, loops[i].seam.x(), loops[i].seam.y(), expected[i].seam.x(), expected[i].seam.y());
        CHECK((loops[i].seam - expected[i].seam).norm() < 0.01);
    }
}

TEST_CASE("The whole-perimeter Blocked warning names the first modifier and counts them all", "[PreciseSeam]")
{
    ExportFixture fixture("back");
    // Separate layers: where both covered the same perimeter, only the higher-priority one would count.
    fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, make_cube(30, 30, 1.4), Vec3d(-5, -5, -1), "first blocker");
    fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, make_cube(30, 30, 2), Vec3d(-5, -5, 1), "second blocker");
    fixture.export_gcode();
    const std::string warning = fixture.warning();
    CAPTURE(warning);
    CHECK(warning.find("Seam Blocked modifier \"first blocker\" of") != std::string::npos);
    CHECK(warning.find("(2 in total) covers an entire perimeter") != std::string::npos);
    CHECK(warning.find("second blocker") == std::string::npos);
}

TEST_CASE("Belt printers slice Precise Seam modifiers in the frame the object was sliced in", "[PreciseSeam][belt]")
{
    Model model;
    Print print;
    Test::init_print({Test::cube(20)}, print, model, {
        { "belt_printer",               1 },
        { "belt_slice_rotation",        "x" },
        { "belt_slice_rotation_angle",  45 },
        { "layer_height",               0.2 },
        { "initial_layer_print_height", 0.2 },
        { "skirt_loops",                0 },
        // Keep the first layer's islands the raw slice, like the modifier's.
        { "elefant_foot_compensation",  0 },
    });
    print.process();
    const PrintObject *object = print.get_object(0);
    REQUIRE(object->layer_count() > 0);

    // A modifier with the object's own mesh and placement must slice to the object's own islands
    // on every layer. Sliced without the belt rotation it would give 20 mm squares on the lower
    // layers and nothing above 20 mm, while the tilted cube reaches about 28 mm.
    Model modifiers;
    ModelVolume *modifier = modifiers.add_object()->add_volume(*model.objects.front()->volumes.front());
    modifier->set_type(ModelVolumeType::PRECISE_SEAM_BLOCKED);
    const std::vector<ExPolygons> slices = object->slice_single_volume_regions(modifier);
    REQUIRE(slices.size() == object->layer_count());

    const double tolerance = scale_(0.05) * scale_(20.);
    for (size_t i = 0; i < slices.size(); ++i) {
        CAPTURE(i, object->get_layer(int(i))->slice_z);
        const ExPolygons &islands = object->get_layer(int(i))->lslices;
        double object_area = 0., modifier_area = 0.;
        for (const ExPolygon &island : islands)
            object_area += island.area();
        for (const ExPolygon &region : slices[i])
            modifier_area += region.area();
        CHECK(std::abs(modifier_area - object_area) < tolerance);
        CHECK(get_extents(slices[i]).inflated(scale_(0.05)).contains(get_extents(islands)));
    }
}

namespace {
// Seam x of a loop lies on the front-face helper (x = 4..6 mm, so -6..-4 relative to the centre).
bool seam_on_front_helper(const OuterLoop &loop)
{
    return loop.seam.x() >= -6.01 && loop.seam.x() <= -3.99 &&
           std::abs(loop.seam.y() - (loop.bounds.min.y() - loop.bounds.center().y())) < 0.01;
}
} // namespace

TEST_CASE("Layer rules limit a modifier to the layers numbered as on the preview slider", "[PreciseSeam][Rules]")
{
    // With a raft, its two layers take the slider numbers 1 and 2, so the cube starts at layer 3.
    const auto [raft_layers, rule] = GENERATE(table<int, std::string>({{0, "1-3"}, {2, "3-5"}}));
    CAPTURE(raft_layers, rule);
    ExportFixture fixture("back", raft_layers);
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front {ps=" + rule + "=ps}");
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i, loops[i].seam.x(), loops[i].seam.y());
        CHECK(seam_on_front_helper(loops[i]) == (i < 3));
    }
    CHECK(fixture.warning().empty());
}

TEST_CASE("The rules field combines with the name and can retype the modifier", "[PreciseSeam][Rules]")
{
    // A Blocked volume relabelled Center by its field acts as a strong modifier on even layers only.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front {ps=center=ps}");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "excl, odd";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i);
        CHECK(seam_on_front_helper(loops[i]) == (i % 2 == 1)); // Loop i is layer i + 1.
    }
    // The type differs from the object list, which is reported once with the modifier name.
    CHECK(fixture.warning().find("front {ps=center=ps}") != std::string::npos);
}

TEST_CASE("Changing only the rules updates the exported seam", "[PreciseSeam][Rules]")
{
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front");
    const auto all_layers = outer_wall_loops(fixture.export_gcode());
    REQUIRE(all_layers.size() == 10);
    CHECK(std::all_of(all_layers.begin(), all_layers.end(), seam_on_front_helper));

    // The geometry is unchanged, so only the export step can pick up the new rules.
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "off";
    const auto switched_off = outer_wall_loops(fixture.export_gcode());
    REQUIRE(switched_off.size() == 10);
    CHECK(std::none_of(switched_off.begin(), switched_off.end(), seam_on_front_helper));
    // A volume switched off by its rules is not reported as having no effect.
    CHECK(fixture.warning().empty());
    // A blank one acts on no seam either, whatever its type, and is not reported.
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "blank,center";
    const auto blank = outer_wall_loops(fixture.export_gcode());
    REQUIRE(blank.size() == 10);
    CHECK(std::none_of(blank.begin(), blank.end(), seam_on_front_helper));
    CHECK(fixture.warning().empty());

    // A rule problem is reported; removing the last rules returns the object to the path without rules,
    // and the old rule warning goes with them.
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "bad";
    const auto invalid = outer_wall_loops(fixture.export_gcode());
    REQUIRE(invalid.size() == 10);
    CHECK(std::all_of(invalid.begin(), invalid.end(), seam_on_front_helper));
    CHECK_FALSE(fixture.warning().empty());
    fixture.model.objects.front()->volumes.back()->precise_seam_rules.clear();
    const auto restored = outer_wall_loops(fixture.export_gcode());
    REQUIRE(restored.size() == 10);
    CHECK(std::all_of(restored.begin(), restored.end(), seam_on_front_helper));
    CHECK(fixture.warning().empty());
}

TEST_CASE("Without rules the export collects no layer plan early", "[PreciseSeam][Rules]")
{
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front");
    fixture.print.apply(fixture.model, fixture.config);
    // Counts every use of what the export lends; without rules there is none.
    int calls = 0;
    const PreciseSeam::ExportPlan::Access access{
        [&calls] { ++calls; return PreciseSeam::ExportPlan::LayersByZ{}; },
        [&calls](const PrintObject &) { ++calls; return PreciseSeam::ExportPlan::ObjectLayers{}; },
        [&calls](const PrintObject &, const Layer &) { ++calls; return true; }};
    for (const bool by_object : {false, true}) {
        CAPTURE(by_object);
        PreciseSeam::ExportPlan plan;
        plan.prepare(fixture.print, ToolOrdering{}, by_object, access);
        CHECK(calls == 0);
        CHECK(plan.numbering().empty());
        CHECK_FALSE(plan.take_layers_by_z().has_value());
        CHECK(plan.object_layers(*fixture.print.objects().front()) == nullptr);
    }

    // The control: rules make the export collect the plan of the plate once.
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "odd";
    fixture.print.apply(fixture.model, fixture.config);
    PreciseSeam::ExportPlan plan;
    plan.prepare(fixture.print, ToolOrdering{}, false, access);
    CHECK(calls == 1);
    CHECK(plan.take_layers_by_z().has_value());
}

TEST_CASE("Printing by object, an object with rules and one without keep their own seams", "[PreciseSeam][Rules]")
{
    // Only the object with rules has its layer plan collected before the seams; the other one is
    // collected where the export always does it.
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("print_sequence", "by object");
    config.set_deserialize_strict("seam_position", "back");
    config.set_deserialize_strict("seam_slope_type", "none");
    config.set_deserialize_strict("layer_height", "0.2");
    config.set_deserialize_strict("initial_layer_print_height", "0.2");
    config.set_deserialize_strict("outer_wall_line_width", "0.4");
    config.set_deserialize_strict("initial_layer_line_width", "0.4");
    Model model;
    Print print;
    Test::init_print({make_cube(20, 20, 2), make_cube(20, 20, 2)}, print, model, config);
    REQUIRE(model.objects.size() == 2);
    // The same front-face helper as ExportFixture::add() on each object.
    for (ModelObject *object : model.objects) {
        const ModelVolume *cube = object->volumes.front();
        const Vec3d cube_min = cube->mesh().transformed_bounding_box(cube->get_matrix()).min;
        ModelVolume *helper = object->add_volume(make_cube(2, 6, 4));
        helper->set_type(ModelVolumeType::PRECISE_SEAM_CENTER);
        helper->name = "front";
        const Vec3d current = helper->mesh().transformed_bounding_box(helper->get_matrix()).min;
        helper->set_offset(helper->get_offset() + cube_min + Vec3d(4, -3, -1) - current);
    }
    model.objects.front()->volumes.back()->precise_seam_rules = "odd";
    print.apply(model, config);
    const auto loops = outer_wall_loops(Test::gcode(print));
    REQUIRE(loops.size() == 20);
    // The objects print one after the other, ten loops each, in an order the arrangement decides. Each one
    // numbers its own layers from 1: the object with rules on its odd layers, the other on every layer.
    const auto on_helper = [&loops](size_t first) {
        std::vector<bool> out;
        for (size_t i = first; i < first + 10; ++i)
            out.push_back(seam_on_front_helper(loops[i]));
        return out;
    };
    std::vector<bool> odd_layers;
    for (size_t i = 0; i < 10; ++i)
        odd_layers.push_back(i % 2 == 0); // Loop i is layer i + 1.
    const std::vector<bool> first = on_helper(0), second = on_helper(10);
    const bool rules_first = first == odd_layers;
    CHECK((rules_first ? first : second) == odd_layers);
    CHECK((rules_first ? second : first) == std::vector<bool>(10, true));
    CHECK(print.step_state_with_warnings(psGCodeExport).warnings.empty());
}

TEST_CASE("Rules that keep every layer and type give the seams of no rules", "[PreciseSeam][Rules]")
{
    // Objects with rules take their own branch; rules that change nothing must match the branch without them.
    const auto [strong_rules, weak_rules] = GENERATE(table<std::string, std::string>({{"all", "all"}, {"center, all", "green"}}));
    CAPTURE(strong_rules, weak_rules);
    // A strong helper on the front face of the lower layers, an enforcer on the right face of all layers.
    const auto export_loops = [](const std::string &strong, const std::string &weak) {
        ExportFixture fixture("back");
        fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 2), Vec3d(4, -3, -1), "front");
        fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, make_cube(6, 2, 4), Vec3d(17, 4, -1), "right");
        const ModelVolumePtrs &volumes = fixture.model.objects.front()->volumes;
        volumes[1]->precise_seam_rules = strong;
        volumes[2]->precise_seam_rules = weak;
        std::vector<OuterLoop> loops = outer_wall_loops(fixture.export_gcode());
        CHECK(fixture.warning().empty());
        return loops;
    };
    const std::vector<OuterLoop> without_rules = export_loops({}, {});
    REQUIRE(without_rules.size() == 10);
    // The strong helper acts on the lower layers only, so both kinds of helper take part.
    CHECK(seam_on_front_helper(without_rules.front()));
    CHECK_FALSE(seam_on_front_helper(without_rules.back()));

    const std::vector<OuterLoop> with_rules = export_loops(strong_rules, weak_rules);
    REQUIRE(with_rules.size() == without_rules.size());
    for (size_t i = 0; i < with_rules.size(); ++i) {
        CAPTURE(i);
        CHECK(with_rules[i].seam.x() == without_rules[i].seam.x());
        CHECK(with_rules[i].seam.y() == without_rules[i].seam.y());
    }
}

TEST_CASE("Each layer rule problem is a separate export warning", "[PreciseSeam][Rules]")
{
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front {ps=bad, 1-3=ps}");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "500";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    // "500" selects no layer of ten and is applied: the modifier then acts on layers 1-3 only.
    for (size_t i = 0; i < loops.size(); ++i)
        CHECK(seam_on_front_helper(loops[i]) == (i < 3));
    std::vector<std::string> rule_warnings;
    for (const auto &warning : fixture.print.step_state_with_warnings(psGCodeExport).warnings)
        if (warning.message.find("rule \"") != std::string::npos)
            rule_warnings.push_back(warning.message);
    REQUIRE(rule_warnings.size() == 2);
    CHECK(rule_warnings[0].find("\"bad\"") != std::string::npos);
    CHECK(rule_warnings[1].find("\"500\"") != std::string::npos);
    for (const std::string &warning : rule_warnings)
        CHECK(warning.find('\n') == std::string::npos);
}

namespace {
// Rule warnings of the last export, in issue order.
std::vector<std::string> rule_warnings(const ExportFixture &fixture)
{
    std::vector<std::string> out;
    for (const auto &warning : fixture.print.step_state_with_warnings(psGCodeExport).warnings)
        // The combined Precise Seam warning has its own id; rule warnings use the default one.
        if (warning.message_id == 0 && warning.message.find("Precise Seam: modifier") != std::string::npos)
            out.push_back(warning.message);
    return out;
}
} // namespace

// Hidden: the base export crashes in this configuration as well (SIGSEGV at 501139e527, without Precise
// Seam), so the test cannot run until that is fixed. It guards the plan choice of the layer rules.
TEST_CASE("A by-object print exported by layer for its wipe tower keeps its layers and rules", "[PreciseSeam][Rules][Regression][.]")
{
    // One object printed by object with a smooth timelapse gets a wipe tower, and the export then
    // runs by layer: the layer plan and its numbers must follow that branch.
    ExportFixture fixture("back");
    fixture.config.set_deserialize_strict("print_sequence", "by object");
    fixture.config.set_deserialize_strict("enable_prime_tower", "1");
    fixture.config.set_deserialize_strict("timelapse_type", "1");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front {ps=1-3=ps}");
    const std::string gcode = fixture.export_gcode();
    // The tower is printed, so the export took the by-layer branch.
    REQUIRE(gcode.find("TYPE:Prime tower") != std::string::npos);
    const auto loops = outer_wall_loops(gcode);
    REQUIRE(loops.size() == 10);
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i);
        CHECK(seam_on_front_helper(loops[i]) == (i < 3));
    }
}

TEST_CASE("Rule warnings stay distinct per source and count equal texts", "[PreciseSeam][Rules]")
{
    ExportFixture fixture("back");
    // Two modifiers with the same name and the same bad field give one message with its count.
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, make_cube(2, 6, 4), Vec3d(4, -3, -1), "helper");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "bad";
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, make_cube(2, 6, 4), Vec3d(14, -3, -1), "helper");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "bad";
    // One modifier with the same bad rule in its name and its field gives two messages.
    fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, make_cube(2, 6, 4), Vec3d(9, -3, -1), "twice {ps=bad=ps}");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "bad";
    fixture.export_gcode();
    const std::vector<std::string> warnings = rule_warnings(fixture);
    REQUIRE(warnings.size() == 3);
    CHECK(warnings[0].find("\"helper\"") != std::string::npos);
    CHECK(warnings[0].find("(2 in total)") != std::string::npos);
    CHECK(warnings[1].find("name block 1, element 1") != std::string::npos);
    CHECK(warnings[2].find("rules field, element 1") != std::string::npos);
}

TEST_CASE("Rule warnings keep one line when the rule text holds a line break", "[PreciseSeam][Rules]")
{
    // "50<LF>0" is the valid set 500, beyond the ten layers; the message must keep its reason.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "50\n0";
    fixture.export_gcode();
    const std::vector<std::string> warnings = rule_warnings(fixture);
    REQUIRE(warnings.size() == 1);
    CHECK(warnings.front().find('\n') == std::string::npos);
    CHECK(warnings.front().find('\r') == std::string::npos);
    CHECK(warnings.front().find("50\\n0") != std::string::npos);
    CHECK(warnings.front().find("selects no existing layer") != std::string::npos);
}

TEST_CASE("Rule warnings of one modifier follow the element order", "[PreciseSeam][Rules]")
{
    // The empty set of element 1 is found after parsing, yet it is reported before element 2.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "500,bad";
    fixture.export_gcode();
    const std::vector<std::string> warnings = rule_warnings(fixture);
    REQUIRE(warnings.size() == 2);
    CHECK(warnings[0].find("\"500\"") != std::string::npos);
    CHECK(warnings[1].find("\"bad\"") != std::string::npos);
}

TEST_CASE("Height rules keep their place in mm when the layer height changes", "[PreciseSeam][Rules][Heights]")
{
    // Loop i is the layer (i h, (i + 1) h]: z0.5-1.2 takes the layers that meet the range, zin0.5-1.2 those inside it.
    const auto [height, rule, first, last] = GENERATE(table<std::string, std::string, size_t, size_t>(
        {{"0.2", "z0.5-1.2", 2, 5}, {"0.1", "z0.5-1.2", 4, 11}, {"0.2", "zin0.5-1.2", 3, 5}, {"0.1", "zin0.5-1.2", 5, 11}}));
    CAPTURE(height, rule);
    ExportFixture fixture("back");
    fixture.config.set_deserialize_strict("layer_height", height);
    fixture.config.set_deserialize_strict("initial_layer_print_height", height);
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front {ps=" + rule + "=ps}");
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == (height == "0.2" ? 10 : 20));
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i);
        CHECK(seam_on_front_helper(loops[i]) == (i >= first && i <= last));
    }
    CHECK(fixture.warning().empty());
}

TEST_CASE("Heights count from the bed, so a raft lifts the model layers", "[PreciseSeam][Rules][Heights]")
{
    ExportFixture fixture("back", 2);
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front {ps=z0-2=ps}");
    const auto loops = outer_wall_loops(fixture.export_gcode());
    const PrintObject &object = *fixture.print.objects().front();
    REQUIRE(loops.size() == 10);
    REQUIRE(object.layers().size() == loops.size());
    size_t selected = 0;
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i, object.layers()[i]->bottom_z());
        const bool expected = object.layers()[i]->bottom_z() < 2. - EPSILON;
        CHECK(seam_on_front_helper(loops[i]) == expected);
        selected += expected;
    }
    // Without the raft every layer of the 2 mm cube starts below 2 mm.
    CHECK(selected > 0);
    CHECK(selected < loops.size());
    CHECK(fixture.warning().empty());
}

TEST_CASE("override in the field replaces the name rules of the exported seam", "[PreciseSeam][Rules][Override]")
{
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front {ps=1-3=ps}");
    ModelVolume &helper = *fixture.model.objects.front()->volumes.back();
    helper.precise_seam_rules = "override, excl, odd";
    const auto replaced = outer_wall_loops(fixture.export_gcode());
    REQUIRE(replaced.size() == 10);
    for (size_t i = 0; i < replaced.size(); ++i) {
        CAPTURE(i);
        CHECK(seam_on_front_helper(replaced[i]) == (i % 2 == 1)); // Loop i is layer i + 1.
    }
    CHECK(fixture.warning().empty());

    // Clearing the field brings the name rules back.
    helper.precise_seam_rules.clear();
    const auto restored = outer_wall_loops(fixture.export_gcode());
    REQUIRE(restored.size() == 10);
    for (size_t i = 0; i < restored.size(); ++i) {
        CAPTURE(i);
        CHECK(seam_on_front_helper(restored[i]) == (i < 3));
    }
}

namespace {
// Index of the front-face helper whose centre x (relative to the loop centre) holds the seam, or -1.
int front_helper(const OuterLoop &loop, const std::vector<double> &centres)
{
    if (std::abs(loop.seam.y() - (loop.bounds.min.y() - loop.bounds.center().y())) >= 0.01)
        return -1;
    for (size_t i = 0; i < centres.size(); ++i)
        if (std::abs(loop.seam.x() - centres[i]) < 0.01)
            return int(i);
    return -1;
}

std::vector<int> front_helpers(const std::vector<OuterLoop> &loops, const std::vector<double> &centres)
{
    std::vector<int> out;
    for (const OuterLoop &loop : loops)
        out.push_back(front_helper(loop, centres));
    return out;
}

// Front-face helper 2 x 6 mm that ends `top` mm above the bottom of the object; placed 1 mm below it.
TriangleMesh front_box(double top = 3.) { return make_cube(2, 6, top + 1.); }

// Two 8 x 20 x 2 mm islands of one object, L at x = 0..8 and R at x = 12..20; each layer has both loops.
TriangleMesh two_islands()
{
    TriangleMesh mesh  = make_cube(8, 20, 2);
    TriangleMesh right = make_cube(8, 20, 2);
    right.translate(12., 0., 0.);
    mesh.merge(right);
    return mesh;
}

// Seam helpers of the L and R loops per layer, the left island first.
std::vector<std::pair<int, int>> island_helpers(const std::vector<OuterLoop> &loops, const std::vector<double> &centres)
{
    std::vector<std::pair<int, int>> out;
    for (size_t i = 0; i + 1 < loops.size(); i += 2) {
        const bool       left_first = loops[i].bounds.center().x() < loops[i + 1].bounds.center().x();
        const OuterLoop &left = loops[left_first ? i : i + 1], &right = loops[left_first ? i + 1 : i];
        out.emplace_back(front_helper(left, centres), front_helper(right, centres));
    }
    return out;
}
} // namespace

TEST_CASE("Relay members take the seam in turn, layer by layer", "[PreciseSeam][Rules][Relay]")
{
    // A and C are Center at x = 3 and 17 mm, B is Left at x = 9 mm: -7, -1 and 7 relative to the centre.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_LEFT, front_box(), Vec3d(9, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(16, -3, -1), "C {ps=relay=ps}");
    const std::vector<double> centres{-7., -1., 7.};
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_helpers(loops, centres) == std::vector<int>{0, 1, 2, 0, 1, 2, 0, 1, 2, 0});
    CHECK(fixture.warning().empty());
    // The schedule is rebuilt for every export and comes out the same.
    CHECK(front_helpers(outer_wall_loops(fixture.export_gcode()), centres) == front_helpers(loops, centres));
}

TEST_CASE("A relay member that cannot leaves the turn with it", "[PreciseSeam][Rules][Relay]")
{
    // B reaches layers 1-5 only by geometry, or acts on layers 6-10 only by its rules.
    const auto [top, rules, expected] = GENERATE(table<double, std::string, std::vector<int>>({
        {1., "relay", {0, 1, 0, 1, 0, 0, 0, 0, 0, 0}},
        {3., "relay,6-10", {0, 0, 0, 0, 0, 1, 0, 1, 0, 1}}}));
    CAPTURE(rules);
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(top), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = rules;
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_helpers(loops, {-7., 7.}) == expected);
    CHECK(fixture.warning().empty());
}

TEST_CASE("A senior modifier takes the layer without spending the relay turn", "[PreciseSeam][Rules][Relay]")
{
    // X, first in the object list and without relay, reaches layer 1 only; A then takes the first turn.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(0.2), Vec3d(9, -3, -1), "X");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_helpers(loops, {0., -7., 7.}) == std::vector<int>{0, 1, 2, 1, 2, 1, 2, 1, 2, 1});
}

TEST_CASE("A relay group acts at the place of its most senior member", "[PreciseSeam][Rules][Relay]")
{
    // A, X, B in the object list: the block of A and B stands before X, so X never takes the seam.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay2";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(9, -3, -1), "X");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay2";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_helpers(loops, {-7., 0., 7.}) == std::vector<int>{0, 2, 0, 2, 0, 2, 0, 2, 0, 2});
    // X was never evaluated, so it is not reported as having no effect.
    CHECK(fixture.warning().empty());
}

TEST_CASE("Another relay member stands in where the holder cannot, without taking the turn", "[PreciseSeam][Rules][Relay]")
{
    // M crosses the front of both islands (x = 3 and 15 mm), N only the left one (x = 6 mm). Relative to the
    // island centres (x = 4 and 16 mm) these are -1 and 2.
    ExportFixture fixture("back", 0, two_islands());
    TriangleMesh  both = front_box(), right = front_box();
    right.translate(12., 0., 0.);
    both.merge(right);
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, std::move(both), Vec3d(2, -3, -1), "M");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(5, -3, -1), "N");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 20);
    // On even layers N holds the turn on the left island, and M stands in on the right one.
    std::vector<std::pair<int, int>> expected;
    for (int layer = 1; layer <= 10; ++layer)
        expected.emplace_back(layer % 2 == 1 ? 0 : 1, 0);
    CHECK(island_helpers(loops, {-1., 2.}) == expected);
    CHECK(fixture.warning().empty());
}

TEST_CASE("Relay members on different islands each place their seam on every layer", "[PreciseSeam][Rules][Relay]")
{
    // A crosses only the left island (x = 3 mm), B only the right one (x = 18 mm): -1 and 2 relative to the
    // island centres. One group covers both, as each member stands in for the other.
    ExportFixture fixture("back", 0, two_islands());
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(17, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 20);
    CHECK(island_helpers(loops, {-1., 2.}) == std::vector<std::pair<int, int>>(10, {0, 1}));
}

TEST_CASE("A relay group that mixes strong and weak members acts without relay and is reported", "[PreciseSeam][Rules][Relay]")
{
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1,red";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    // A acts as an ordinary strong modifier on every layer.
    CHECK(front_helpers(loops, {-7., 7.}) == std::vector<int>(10, 0));
    const std::vector<std::string> warnings = rule_warnings(fixture);
    REQUIRE(warnings.size() == 3); // B also acts with the type from its rules.
    CHECK(warnings[0].find("\"A\"") != std::string::npos);
    CHECK(warnings[0].find("\"relay1\" acts without relay") != std::string::npos);
    CHECK(warnings[1].find("\"B\"") != std::string::npos);
    CHECK(warnings[2].find("\"B\"") != std::string::npos);
    CHECK(warnings[2].find("\"relay1\" acts without relay") != std::string::npos);
}

TEST_CASE("A stand-in tries the members after the holder first, round the queue", "[PreciseSeam][Rules][Relay]")
{
    // A and C cross only the left island (x = 3 and 6 mm, -1 and 2 relative to its centre), B only the right
    // one (x = 15 mm, -1). While B holds, the left island goes to C, the member after B, not to A.
    ExportFixture fixture("back", 0, two_islands());
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(14, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(5, -3, -1), "C");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 20);
    // Holders A, B, C, A, ...: the left island shows A, then C under B and under C. Indices into {A, C}.
    std::vector<std::pair<int, int>> expected;
    for (int layer = 1; layer <= 10; ++layer)
        expected.emplace_back(layer % 3 == 1 ? 0 : 1, 0);
    CHECK(island_helpers(loops, {-1., 2.}) == expected);
}

TEST_CASE("A stand-in of a senior group keeps the turn of the junior group", "[PreciseSeam][Rules][Relay]")
{
    // Object list H_A, G_C, H_B, G_D: groups H = relay1 (anchor H_A) and G = relay2 (anchor G_C) interleave.
    // H_A crosses only the left island; H_B, G_C and G_D only the right one, H_B on odd layers by its rules.
    // Relative to the right island's centre (x = 16 mm) H_B, G_C and G_D are at -2.5, 0 and 2.5.
    ExportFixture fixture("back", 0, two_islands());
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "H_A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(15, -3, -1), "G_C");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay2";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(12.5, -3, -1), "H_B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1,odd";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(17.5, -3, -1), "G_D");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay2";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 20);
    // H takes the right island on odd layers, as a stand-in on layer 1 and as the holder on layer 3, so G
    // has no holder there and keeps its turn: G_C on layer 2, G_D on layer 4. Were the stand-in left out of
    // the seniors of G, G would spend its turn on layer 1 and show G_D on layer 2.
    // Indices into {H_A, H_B, G_C, G_D}.
    const std::vector<std::pair<int, int>> expected{{0, 1}, {0, 2}, {0, 1}, {0, 3}, {0, 1},
                                                    {0, 2}, {0, 1}, {0, 3}, {0, 1}, {0, 2}};
    CHECK(island_helpers(loops, {-1., -2.5, 0., 2.5}) == expected);
}

// Expected holders of relay-random come from an independent model of the hash in SPEC v7.
TEST_CASE("relay-random members take the seam at random, never twice in a row", "[PreciseSeam][Rules][Relay]")
{
    // A, B and C as in the relay test above: -7, -1 and 7 relative to the centre; each wins on every layer.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random";
    fixture.add(ModelVolumeType::PRECISE_SEAM_LEFT, front_box(), Vec3d(9, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(16, -3, -1), "C {ps=relay-random=ps}");
    const std::vector<double> centres{-7., -1., 7.};
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_helpers(loops, centres) == std::vector<int>{0, 1, 2, 1, 2, 0, 2, 1, 2, 0});
    CHECK(fixture.warning().empty());
    // Nothing is kept between exports: the next one gives the same seams.
    CHECK(front_helpers(outer_wall_loops(fixture.export_gcode()), centres) == front_helpers(loops, centres));
}

TEST_CASE("relay1 and relay-random1 are different groups", "[PreciseSeam][Rules][Relay]")
{
    // H = relay1 crosses only the left island (H_A, H_B at -1 and 2 relative to its centre), G = relay-random1
    // only the right one (G_C, G_D, G_E at -2.5, 0 and 2.5). H, the senior group, leaves the right island
    // free, so both groups show their own order.
    ExportFixture fixture("back", 0, two_islands());
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "H_A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(12.5, -3, -1), "G_C");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(5, -3, -1), "H_B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(15, -3, -1), "G_D");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(17.5, -3, -1), "G_E");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random1";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 20);
    // Indices into {H_A, H_B, G_C, G_D, G_E}: H alternates, G follows its random holders.
    const std::vector<int> g{4, 2, 4, 3, 4, 3, 2, 3, 4, 2};
    std::vector<std::pair<int, int>> expected;
    for (int layer = 1; layer <= 10; ++layer)
        expected.emplace_back(layer % 2 == 1 ? 0 : 1, g[layer - 1]);
    CHECK(island_helpers(loops, {-1., 2., -2.5, 0., 2.5}) == expected);
    CHECK(fixture.warning().empty());
}

TEST_CASE("A layer that a senior modifier takes keeps the previous relay-random holder", "[PreciseSeam][Rules][Relay]")
{
    // X, first in the object list and without relay, acts on layer 3 only (0 relative to the centre); A, B
    // and C at -7, -3 and 7 win on every other layer.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(9, -3, -1), "X");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "3";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(6, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random1";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(16, -3, -1), "C");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random1";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    // Holders C, A, -, C, ...: layer 4 leaves out A, the holder of layer 2; had layer 3 cleared it, B would
    // hold layer 4. Indices into {X, A, B, C}.
    CHECK(front_helpers(loops, {0., -7., -3., 7.}) == std::vector<int>{3, 1, 0, 3, 2, 3, 1, 2, 3, 1});
}

TEST_CASE("A relay-random group that mixes strong and weak members acts without relay and is reported", "[PreciseSeam][Rules][Relay]")
{
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random,red";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_helpers(loops, {-7., 7.}) == std::vector<int>(10, 0));
    const std::vector<std::string> warnings = rule_warnings(fixture);
    REQUIRE(warnings.size() == 3); // B also acts with the type from its rules.
    CHECK(warnings[0].find("\"relay-random\" acts without relay") != std::string::npos);
    CHECK(warnings[2].find("\"relay-random\" acts without relay") != std::string::npos);
}

TEST_CASE("A stand-in of a senior group leaves the relay-random group no holder on its layers", "[PreciseSeam][Rules][Relay]")
{
    // Object list H_A, G_C, H_B, G_D, G_E: H = relay1 (anchor H_A), G = relay-random1 (anchor G_C). H_A crosses
    // only the left island (1.5 relative to its centre); H_B, G_C, G_D and G_E, 1 mm wide, only the right one
    // (-3, -1, 1 and 3), H_B on odd layers by its rules.
    ExportFixture fixture("back", 0, two_islands());
    const auto add = [&fixture](TriangleMesh mesh, double x, const char *name, const char *rules) {
        fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, std::move(mesh), Vec3d(x, -3, -1), name);
        fixture.model.objects.front()->volumes.back()->precise_seam_rules = rules;
    };
    const TriangleMesh narrow = make_cube(1, 6, 4);
    add(front_box(), 4.5, "H_A", "relay1");
    add(narrow, 14.5, "G_C", "relay-random1");
    add(narrow, 12.5, "H_B", "relay1,odd");
    add(narrow, 16.5, "G_D", "relay-random1");
    add(narrow, 18.5, "G_E", "relay-random1");
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 20);
    // H_B takes the right island on odd layers, as a stand-in where H_A holds, so G wins only on even layers
    // and keeps its previous holder over the odd ones: G_D, G_E, G_D, G_C, G_D. Had the seniors of G held
    // only the holder of H, the pre-pass would give layer 1 to G and layer 2 would go to G_C; had an odd layer
    // cleared the previous holder, layer 4 would go to G_D. Indices into {H_A, H_B, G_C, G_D, G_E}.
    const std::vector<int> right{1, 3, 1, 4, 1, 3, 1, 2, 1, 3};
    std::vector<std::pair<int, int>> expected;
    for (int layer = 1; layer <= 10; ++layer)
        expected.emplace_back(0, right[layer - 1]);
    CHECK(island_helpers(loops, {1.5, -3., -1., 1., 3.}) == expected);
}

namespace {
// Index of the front-face zone 2 mm wide from `lefts[i]` mm along the cube that holds the seam, or -1.
int front_zone(const OuterLoop &loop, const std::vector<double> &lefts)
{
    if (std::abs(loop.seam.y() - (loop.bounds.min.y() - loop.bounds.center().y())) >= 0.01)
        return -1;
    for (size_t i = 0; i < lefts.size(); ++i)
        if (loop.seam.x() >= lefts[i] - 10.01 && loop.seam.x() <= lefts[i] - 7.99)
            return int(i);
    return -1;
}

std::vector<int> front_zones(const std::vector<OuterLoop> &loops, const std::vector<double> &lefts)
{
    std::vector<int> out;
    for (const OuterLoop &loop : loops)
        out.push_back(front_zone(loop, lefts));
    return out;
}

// A back seam lies on the rear face (y = 20 mm), 10 mm behind the centre.
bool seam_at_back(const OuterLoop &loop) { return loop.seam.y() > 7.; }
} // namespace

TEST_CASE("Weak relay members act in turn, layer by layer", "[PreciseSeam][Rules][Relay]")
{
    // Enforced zones at x = 2, 9 and 16 mm on the front face; the back seam would lie on the rear face.
    const std::string label = GENERATE(as<std::string>{}, "relay", "relay-random");
    CAPTURE(label);
    ExportFixture fixture("back");
    for (double x : {2., 9., 16.}) {
        fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(x, -3, -1));
        fixture.model.objects.front()->volumes.back()->precise_seam_rules = label;
    }
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    // relay-random takes the holders of the strong relay-random test: every member acts on every layer.
    CHECK(front_zones(loops, {2., 9., 16.}) ==
          (label == "relay" ? std::vector<int>{0, 1, 2, 0, 1, 2, 0, 1, 2, 0} : std::vector<int>{0, 1, 2, 1, 2, 0, 2, 1, 2, 0}));
    CHECK(fixture.warning().empty());
}

TEST_CASE("A weak relay member keeps its turn where its zone ends, and nobody stands in", "[PreciseSeam][Rules][Relay]")
{
    // A at x = 2 mm crosses every layer; B at x = 16 mm ends 1 mm above the bottom, so layers 1-5 only.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(1.), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    // B holds the even layers; above its zone the seam stays at the back, not in A.
    CHECK(front_zones(loops, {2., 16.}) == std::vector<int>{0, 1, 0, 1, 0, -1, 0, -1, 0, -1});
    for (size_t i : {5, 7, 9})
        CHECK(seam_at_back(loops[i]));
}

TEST_CASE("Layers where only one weak relay member acts go to it", "[PreciseSeam][Rules][Relay]")
{
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay,6-10";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_zones(loops, {2., 16.}) == std::vector<int>{0, 0, 0, 0, 0, 1, 0, 1, 0, 1});
}

TEST_CASE("A layer that a strong modifier takes spends the weak relay turn", "[PreciseSeam][Rules][Relay]")
{
    // X, Center at x = 9..11 mm on layer 3 only; A and B, Enforced at x = 2 and 16 mm, in a relay.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(9, -3, -1), "X");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "3";
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    // Layer 3 is the turn of A, so layer 4 goes to B, as without X. Indices into {A, X, B}.
    CHECK(front_zones(loops, {2., 9., 16.}) == std::vector<int>{0, 2, 1, 2, 0, 2, 0, 2, 0, 2});
}

TEST_CASE("A weak relay group acts at the place of its most senior member", "[PreciseSeam][Rules][Relay]")
{
    // Object list A, X, B: A and B Enforced in a relay at x = 2 and 16 mm, X Blocked over the zone of B. The block
    // {A, B} stands at A, above X, so B wins over X on its layers; at its own place below X it would not.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_BLOCKED, make_cube(6, 6, 4), Vec3d(14, -3, -1), "X");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(16, -3, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_zones(loops, {2., 16.}) == std::vector<int>{0, 1, 0, 1, 0, 1, 0, 1, 0, 1});
}

TEST_CASE("A weak relay group may mix Enforced, Blocked and Neutral members", "[PreciseSeam][Rules][Relay]")
{
    // A, Enforced at x = 2 mm, and B over everything from y = 17 mm back, in a relay. Blocked B moves the seam
    // off the back on its layers; Neutral B, the control, leaves it there.
    const auto type = GENERATE(ModelVolumeType::PRECISE_SEAM_BLOCKED, ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    CAPTURE(int(type));
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(2, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(type, make_cube(30, 6, 4), Vec3d(-5, 17, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i, loops[i].seam.x(), loops[i].seam.y());
        if (i % 2 == 0)
            CHECK(front_zone(loops[i], {2.}) == 0);
        else
            CHECK(seam_at_back(loops[i]) == (type == ModelVolumeType::PRECISE_SEAM_NEUTRAL));
    }
    CHECK(fixture.warning().empty());
}

TEST_CASE("A Neutral weak relay holder clears a junior zone on its layers", "[PreciseSeam][Rules][Relay]")
{
    // Object list A, G, X: A Enforced and G Neutral in a relay, X Enforced without one; A and X at x = 9 mm, G
    // over both. The block {A, G} applies after X, so G clears X on its layers and the seam goes back.
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(9, -3, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_NEUTRAL, make_cube(4, 6, 4), Vec3d(8, -3, -1), "G");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(9, -3, -1), "X");
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (size_t i = 0; i < loops.size(); ++i) {
        CAPTURE(i, loops[i].seam.x(), loops[i].seam.y());
        if (i % 2 == 0)
            CHECK(front_zone(loops[i], {9.}) == 0);
        else
            CHECK(seam_at_back(loops[i]));
    }
}

namespace {
using RuleCode = PreciseSeam::RuleDiagnosticCode;

// Weak volumes of the fixture object initialized on `count` layers numbered from 1, in application order.
// No rule here selects by height, so the layers of the print object are not read.
struct WeakInit {
    std::vector<const ModelVolume*>              weak;
    PreciseSeam::EffectiveModifiers              effective;
    std::optional<PreciseSeam::RelayOrder>       relay;
    std::vector<PreciseSeam::ModifierDiagnostic> diagnostics;

    WeakInit(const ExportFixture &fixture, size_t count)
    {
        PreciseSeam::LayerNumbering numbering;
        numbering.count = count;
        for (size_t i = 1; i <= count; ++i)
            numbering.object_layers.push_back(i);
        std::vector<const ModelVolume*> strong;
        bool                            has_strong = false;
        PreciseSeam::RuleContext context = PreciseSeam::prepare_rule_context(strong, weak, has_strong, *fixture.model.objects.front(),
                                                                             &numbering, fixture.print.objects().front()->layers(),
                                                                             diagnostics);
        effective = std::move(context.effective);
        relay     = std::move(context.relay);
    }

    // Layers of the volume as '1' where it acts; empty when it acts on every layer.
    std::string row(const ModelVolume *volume) const
    {
        std::string out;
        for (bool active : effective.at(volume).layers)
            out += active ? '1' : '0';
        return out;
    }

    size_t reports(RuleCode code) const
    {
        return size_t(std::count_if(diagnostics.begin(), diagnostics.end(),
                                    [code](const PreciseSeam::ModifierDiagnostic &d) { return d.diagnostic.code == code; }));
    }
};

// A Precise Seam volume of the fixture object with the given rules field; its place does not matter here.
const ModelVolume *add_modifier(ExportFixture &fixture, ModelVolumeType type, const std::string &rules)
{
    ModelVolume *volume = fixture.model.objects.front()->add_volume(make_cube(2, 2, 2));
    volume->set_type(type);
    volume->precise_seam_rules = rules;
    return volume;
}
constexpr ModelVolumeType green = ModelVolumeType::PRECISE_SEAM_ENFORCED;
} // namespace

TEST_CASE("A weak relay group stands at its most senior member, and each member keeps the layers it holds", "[PreciseSeam][Rules][Relay]")
{
    using Volumes = std::vector<const ModelVolume*>;
    ExportFixture fixture("back");
    SECTION("pr1 moves the anchor, the queue stays in the object list order")
    {
        const ModelVolume *X = add_modifier(fixture, ModelVolumeType::PRECISE_SEAM_BLOCKED, "");
        const ModelVolume *A = add_modifier(fixture, green, "relay");
        const ModelVolume *B = add_modifier(fixture, green, "relay,pr1");
        const WeakInit init(fixture, 4);
        // Priority order B, X, A becomes the block {A, B} and X, so the block is senior to X only by pr1;
        // weak ones apply in reverse.
        CHECK(init.weak == Volumes{X, B, A});
        CHECK(init.row(A) == "1010");
        CHECK(init.row(B) == "0101");
        CHECK(init.row(X).empty());
        CHECK_FALSE(init.relay.has_value()); // The schedule of weak groups needs no pre-pass.
    }
    SECTION("A heavier priority moves the anchor further: pr10 before pr1")
    {
        const ModelVolume *X = add_modifier(fixture, green, "pr1");
        const ModelVolume *A = add_modifier(fixture, green, "relay");
        const ModelVolume *B = add_modifier(fixture, green, "relay,pr10");
        const WeakInit init(fixture, 4);
        // Priority order B, X, A becomes the block {A, B} and X; weak ones apply in reverse.
        CHECK(init.weak == Volumes{X, B, A});
        CHECK(init.row(A) == "1010");
        CHECK(init.row(B) == "0101");
    }
    SECTION("pr-neg puts its member after the volumes without priority; its group stands at another member")
    {
        const ModelVolume *Y = add_modifier(fixture, green, "relay,pr-neg1");
        const ModelVolume *W = add_modifier(fixture, ModelVolumeType::PRECISE_SEAM_BLOCKED, "");
        const ModelVolume *Z = add_modifier(fixture, green, "relay");
        const WeakInit init(fixture, 4);
        // Priority order W, Z, Y: the block {Y, Z} stands at Z, after W; the queue keeps the object list order.
        CHECK(init.weak == Volumes{Z, Y, W});
        CHECK(init.row(Y) == "1010");
        CHECK(init.row(Z) == "0101");
    }
    SECTION("An anchor without layers keeps the place of the group")
    {
        add_modifier(fixture, green, "relay,pr1,500"); // 500 selects no layer of four.
        const ModelVolume *X = add_modifier(fixture, ModelVolumeType::PRECISE_SEAM_BLOCKED, "");
        const ModelVolume *B = add_modifier(fixture, green, "relay");
        const ModelVolume *C = add_modifier(fixture, green, "relay");
        const WeakInit init(fixture, 4);
        CHECK(init.weak == Volumes{X, C, B});
        CHECK(init.row(B) == "1010");
        CHECK(init.row(C) == "0101");
    }
    SECTION("An off member does not make the group mixed")
    {
        const ModelVolume *A = add_modifier(fixture, green, "relay");
        add_modifier(fixture, ModelVolumeType::PRECISE_SEAM_CENTER, "relay,off");
        const ModelVolume *C = add_modifier(fixture, ModelVolumeType::PRECISE_SEAM_NEUTRAL, "relay");
        const WeakInit init(fixture, 4);
        CHECK(init.weak == Volumes{C, A});
        CHECK(init.row(A) == "1010");
        CHECK(init.row(C) == "0101");
        CHECK(init.reports(RuleCode::RelayGroupMixed) == 0);
    }
    SECTION("A strong member makes the group mixed: every member acts on its own and is reported")
    {
        const ModelVolume *A = add_modifier(fixture, green, "relay");
        add_modifier(fixture, ModelVolumeType::PRECISE_SEAM_CENTER, "relay");
        const ModelVolume *C = add_modifier(fixture, green, "relay");
        const WeakInit init(fixture, 4);
        CHECK(init.weak == Volumes{C, A});
        CHECK(init.row(A).empty());
        CHECK(init.row(C).empty());
        CHECK(init.reports(RuleCode::RelayGroupMixed) == 3);
    }
    SECTION("A single member that takes part acts on its own; none leaves the group out")
    {
        const ModelVolume *A = add_modifier(fixture, green, "relay");
        const ModelVolume *B = add_modifier(fixture, green, "relay,500");
        add_modifier(fixture, green, "relay1,500");
        add_modifier(fixture, green, "relay1,500");
        const WeakInit init(fixture, 4);
        CHECK(init.weak == Volumes{A});
        CHECK(init.row(A).empty());
        CHECK(init.effective.count(B) == 0);
    }
    SECTION("A member that holds no layer is left out")
    {
        const ModelVolume *A = add_modifier(fixture, green, "relay");
        add_modifier(fixture, green, "relay,1");
        add_modifier(fixture, green, "relay,1");
        const WeakInit init(fixture, 4);
        // A holds layer 1 and then is the only member left; the others never get the turn.
        CHECK(init.weak == Volumes{A});
        CHECK(init.row(A) == "1111");
    }
}

// Expected holders come from the independent model of the SPEC v7 hash, over the full table of the active
// members: [2, 0, -, 0, 2, 1, 0, 1, 2, 0]. A table without layer 3 would give other holders from layer 7 on.
TEST_CASE("A weak relay-random group hashes the real layer index and keeps the only active member", "[PreciseSeam][Rules][Relay]")
{
    ExportFixture      fixture("back");
    const ModelVolume *A = add_modifier(fixture, green, "relay-random1,1-2,4-10");
    const ModelVolume *B = add_modifier(fixture, green, "relay-random1,1-2,5-10");
    const ModelVolume *C = add_modifier(fixture, green, "relay-random1,1-2,5-10");
    const WeakInit     init(fixture, 10);
    // No member acts on layer 3; A alone acts on layer 4 and holds it after holding layer 2.
    CHECK(init.row(A) == "0101001001");
    CHECK(init.row(B) == "0000010100");
    CHECK(init.row(C) == "1000100010");
}

namespace {
// Magnet point of the perimeter to the given modifier regions; the limit is in mm.
std::optional<Point> magnet_point(const Polygon &perimeter, ExPolygons regions, std::optional<double> limit = std::nullopt,
                                  bool apex = false)
{
    const std::optional<PreciseSeam::MagnetPoint> found = PreciseSeam::find_magnet_point(
        PreciseSeam::PreparedPerimeter(perimeter), PreciseSeam::prepare_modifier_regions(std::move(regions)),
        limit ? std::optional<double>(scale_(*limit)) : std::nullopt, apex);
    return found ? std::optional<Point>(found->point) : std::nullopt;
}
std::optional<Point> magnet_point(const Polygon &perimeter, const Polygons &loops, std::optional<double> limit = std::nullopt,
                                  bool apex = false)
{
    ExPolygons regions;
    for (const Polygon &loop : loops)
        regions.emplace_back(loop);
    return magnet_point(perimeter, std::move(regions), limit, apex);
}

// The same contour starting at another vertex.
Polygon rotated(Polygon polygon, size_t shift)
{
    std::rotate(polygon.points.begin(), polygon.points.begin() + shift, polygon.points.end());
    return polygon;
}

// An L of 20 x 20 mm whose mass lies on the right: its centroid is at (11.67, 8.33).
Polygon right_heavy_l() { return Polygon(Points{mm(0, 0), mm(20, 0), mm(20, 20), mm(10, 20), mm(10, 10), mm(0, 10)}); }

// The effective rules of a magnet modifier, as prepare_rule_context() sets them.
PreciseSeam::EffectiveModifier magnet_rule(ModelVolumeType type, std::optional<double> limit = std::nullopt)
{
    PreciseSeam::EffectiveModifier rule{type, {}};
    rule.magnet = true;
    if (limit)
        rule.magnet_limit = scale_(*limit);
    return rule;
}
} // namespace

TEST_CASE("A magnet takes the point of the perimeter nearest to the modifier, within its limit", "[PreciseSeam][Magnet]")
{
    const Polygon perimeter = rectangle(0, 0, 20, 20);
    // 3 mm in front of the front face; the parallel edges meet at the projection of the loop centroid.
    const Polygon box = rectangle(8, -5, 12, -3);
    CHECK(magnet_point(perimeter, Polygons{box}) == mm(10, 0));
    CHECK(magnet_point(perimeter, Polygons{box}, 5.) == mm(10, 0));
    CHECK(magnet_point(perimeter, Polygons{box}, 3.) == mm(10, 0));
    // The limit allows the snapping tolerance of 1000 units, 1 um at the default scale, not more.
    CHECK(magnet_point(perimeter, Polygons{box}, 2.9995) == mm(10, 0));
    CHECK_FALSE(magnet_point(perimeter, Polygons{box}, 2.998).has_value());
    CHECK_FALSE(magnet_point(perimeter, Polygons{box}, 2.).has_value());
    // Nested either way, the boundaries decide: a modifier inside the body, and a body inside the modifier.
    CHECK(magnet_point(perimeter, Polygons{rectangle(2, 3, 6, 7)}) == mm(0, 5));
    CHECK(magnet_point(perimeter, Polygons{rectangle(-3, -2, 25, 30)}) == mm(11, 0));
    // A ring around the body: its hole is the nearest boundary.
    ExPolygon ring(rectangle(-10, -10, 30, 30));
    Polygon   hole = rectangle(-2, -3, 23, 25);
    hole.reverse();
    ring.holes.push_back(hole);
    CHECK(magnet_point(perimeter, ExPolygons{ring}) == mm(0, 11));
    CHECK_FALSE(magnet_point(perimeter, ExPolygons{}).has_value());
}

TEST_CASE("magnet0 takes a touch only", "[PreciseSeam][Magnet]")
{
    const Polygon perimeter = rectangle(0, 0, 20, 20);
    // A vertex on the front face, and a vertex on an inclined edge, which computes slightly above zero.
    CHECK(magnet_point(perimeter, Polygons{Polygon(Points{mm(10, 0), mm(8, -2), mm(12, -2)})}, 0.) == mm(10, 0));
    const Polygon triangle(Points{mm(0, 0), mm(20, 0), mm(0, 20)});
    const auto    touch = magnet_point(triangle, Polygons{Polygon(Points{mm(10, 10), mm(14, 12), mm(12, 14)})}, 0.);
    REQUIRE(touch.has_value());
    CHECK((*touch - mm(10, 10)).cast<double>().norm() < 2.);
    // 2 um away is no touch.
    CHECK_FALSE(magnet_point(perimeter, Polygons{Polygon(Points{mm(10, -0.002), mm(8, -2), mm(12, -2)})}, 0.).has_value());
}

TEST_CASE("Equally near magnet points go to the loop centroid, then the perimeter centroid, then rear and left", "[PreciseSeam][Magnet]")
{
    const Polygon perimeter = rectangle(0, 0, 20, 20);
    const size_t  shift     = GENERATE(0, 1, 3);
    CAPTURE(shift);
    SECTION("A vertex at one distance from two edges: its loop centroid picks one, whichever the tree returns")
    {
        // The vertex (2, 2) is 2 mm from the left and the front face; the centroid (4.67, 3.33) is nearer (2, 0).
        const Polygon triangle(Points{mm(2, 2), mm(9, 3), mm(3, 5)});
        CHECK(magnet_point(rotated(perimeter, shift), Polygons{triangle}) == mm(2, 0));
        Polygon mirrored = triangle;
        for (Point &p : mirrored.points)
            p = Point(p.y(), p.x());
        mirrored.reverse();
        CHECK(magnet_point(rotated(perimeter, shift), Polygons{mirrored}) == mm(0, 2));
    }
    SECTION("Each candidate measures to its own loop, so a far region does not move the choice")
    {
        // Equal loops 2 mm in front give (5, 0) and (15, 0); equal for the perimeter centroid too, so the left one.
        const Polygons loops{rectangle(4, -3, 6, -2), rectangle(14, -3, 16, -2)};
        CHECK(magnet_point(rotated(perimeter, shift), loops) == mm(5, 0));
        CHECK(magnet_point(rotated(perimeter, shift), Polygons{loops[0], loops[1], rectangle(40, -30, 60, -20)}) == mm(5, 0));
    }
    SECTION("The perimeter centroid breaks a tie of the loop centroids")
    {
        const Polygons loops{rectangle(2, -3, 4, -2), rectangle(16, -3, 18, -2)};
        CHECK(magnet_point(rotated(right_heavy_l(), shift), loops) == mm(17, 0));
    }
    SECTION("A full symmetry goes to the rear point")
    {
        CHECK(magnet_point(rotated(perimeter, shift), Polygons{rectangle(8, 8, 12, 12)}) == mm(10, 20));
    }
    SECTION("A loop without area has no centroid: the perimeter centroid decides")
    {
        // The flat loop runs 3 mm in front; of the projections of its vertices, x = 9 mm is nearest (10, 10).
        CHECK(magnet_point(rotated(perimeter, shift), Polygons{Polygon(Points{mm(6, -3), mm(12, -3), mm(9, -3)})}) == mm(9, 0));
    }
    SECTION("A candidate within the tolerance of the minimum but beyond the limit is left out")
    {
        // B is 0.8 um farther than A, so both are minima; the perimeter centroid prefers B, the limit only A.
        const Polygons loops{rectangle(2, -3, 4, -2), rectangle(16, -3.0008, 18, -2.0008)};
        CHECK(magnet_point(rotated(right_heavy_l(), shift), loops) == mm(17, 0));
        CHECK(magnet_point(rotated(right_heavy_l(), shift), loops, 1.9995) == mm(3, 0));
    }
}

TEST_CASE("A centroid projection within the limit is found although the tree truncates the query point", "[PreciseSeam][Magnet]")
{
    // Coordinates in units. The perimeter edge lies on 2x + 3y = 4164 and the nearest modifier edge on
    // 2x + 3y = 4164000. The loop centroid (1041, 1041) projects onto the perimeter at (11451, 10410) / 13,
    // 4159836 / sqrt(13) = 1153730.92 units from the modifier, within the limit plus the 1000-unit tolerance.
    // Truncated to (880, 800) it is 1153732.03 away, beyond a margin of one unit.
    Polygon       perimeter(Points{Point(-310218, 208200), Point(0, -100000), Point(314382, -208200)});
    const Polygon loop(Points{Point(-10410000, 8328000), Point(-4160877, 3123), Point(14574000, -8328000)});
    const size_t  shift = GENERATE(0, 1, 2);
    CAPTURE(shift);
    perimeter = rotated(perimeter, shift);
    const auto found = PreciseSeam::find_magnet_point(PreciseSeam::PreparedPerimeter(perimeter),
                                                      PreciseSeam::prepare_modifier_regions(ExPolygons{ExPolygon(loop)}), 1152731.);
    REQUIRE(found.has_value());
    CHECK(found->point == Point(881, 801));
}

TEST_CASE("A magnet with apex moves from the nearest point to the tip of the corner the modifier stands in", "[PreciseSeam][Magnet]")
{
    // 40 x 20 mm with a V notch from the top, its tip at (20, 10). The 1 mm square stands in the notch, near
    // its right side: the nearest point lies on that side, and the tip is 6.2 mm from the square's centre.
    const Polygon notched(Points{mm(0, 0), mm(40, 0), mm(40, 20), mm(25, 20), mm(20, 10), mm(15, 20), mm(0, 20)});
    const Polygon square = rectangle(21, 15.5, 22, 16.5);
    const size_t  shift  = GENERATE(0, 3, 5);
    CAPTURE(shift);
    const Polygon perimeter = rotated(notched, shift);
    const std::optional<Point> plain = magnet_point(perimeter, Polygons{square});
    REQUIRE(plain.has_value());
    CHECK(*plain != mm(20, 10));
    CHECK(magnet_point(perimeter, Polygons{square}, std::nullopt, true) == mm(20, 10));
    // The ray to the far side meets it 3.8 mm from the square: a limit of 5 mm reaches it, 2 mm does not and
    // keeps the nearest point.
    CHECK(magnet_point(perimeter, Polygons{square}, 5., true) == mm(20, 10));
    CHECK(magnet_point(perimeter, Polygons{square}, 2., true) == plain);
    // Inside the body, near the tip of a sharp triangle, the tip wins the same way.
    const Polygon triangle(Points{mm(0, 0), mm(40, 0), mm(20, 40)});
    CHECK(magnet_point(triangle, Polygons{rectangle(20.5, 32.5, 21.5, 33.5)}, std::nullopt, true) == mm(20, 40));
}

TEST_CASE("A magnet with apex keeps the nearest point where no corner holds the modifier", "[PreciseSeam][Magnet]")
{
    const Polygon perimeter = rectangle(0, 0, 20, 20);
    // In front of a face the rays run along it and meet nothing.
    CHECK(magnet_point(perimeter, Polygons{rectangle(8, -5, 12, -3)}, std::nullopt, true) == mm(10, 0));
    // Outside a corner the nearest point is the corner itself.
    CHECK(magnet_point(perimeter, Polygons{rectangle(21, 21, 22, 22)}, std::nullopt, true) == mm(20, 20));
    // Off the centre of a round hole the farthest point of the part is its end, not a vertex inside it.
    Polygon hole;
    for (int i = 0; i < 64; ++i)
        hole.points.push_back(mm(10 + 5 * std::cos(2 * PI * i / 64), 10 + 5 * std::sin(2 * PI * i / 64)));
    const Polygon dot = rectangle(11.4, 9.9, 11.6, 10.1);
    CHECK(magnet_point(hole, Polygons{dot}, std::nullopt, true) == magnet_point(hole, Polygons{dot}));
}

namespace {
// The tip that magnet apex finds between the crossings of the perimeter with the modifier loops.
std::optional<Point> crossing_apex(const Polygon &perimeter, const Polygons &loops)
{
    ExPolygons regions;
    for (const Polygon &loop : loops)
        regions.emplace_back(loop);
    const PreciseSeam::PreparedPerimeter prepared(perimeter);
    const PreciseSeam::ModifierRegions   modifier  = PreciseSeam::prepare_modifier_regions(std::move(regions));
    const PreciseSeam::SegmentExtraction extracted = PreciseSeam::extract_perimeter_segments(prepared, modifier,
                                                                                             ModelVolumeType::PRECISE_SEAM_ENFORCED);
    const std::optional<PreciseSeam::MagnetPoint> found = PreciseSeam::find_crossing_apex(prepared, modifier, extracted);
    return found ? std::optional<Point>(found->point) : std::nullopt;
}

// 40 x 20 mm with a V notch from the top, its tip at (20, 10) and its sides at x = 20 -+ (y - 10) / 2.
Polygon notched() { return Polygon(Points{mm(0, 0), mm(40, 0), mm(40, 20), mm(25, 20), mm(20, 10), mm(15, 20), mm(0, 20)}); }
} // namespace

TEST_CASE("A magnet with apex that crosses both sides of a corner takes its tip", "[PreciseSeam][Magnet]")
{
    const size_t shift = GENERATE(0, 3, 5);
    CAPTURE(shift);
    const Polygon perimeter = rotated(notched(), shift);
    // Across both sides above the tip: the middles of the crossings are (18, 14) and (22, 14), 2 mm from the
    // centroid, and the tip 4 mm; the short part between them holds the tip, the long one the whole body.
    const Polygon across = rectangle(17, 13, 23, 15);
    CHECK(crossing_apex(perimeter, {across}) == mm(20, 10));
    // Over the tip the modifier crosses once, and its type acts.
    CHECK_FALSE(crossing_apex(perimeter, {rectangle(18, 8, 22, 12)}).has_value());
    // So it does when another loop of the modifier crosses once, here the bottom edge.
    CHECK_FALSE(crossing_apex(perimeter, {across, rectangle(5, -1, 7, 1)}).has_value());
    // A loop that does not cross takes no part; the crossing one still finds the tip.
    CHECK(crossing_apex(perimeter, {across, rectangle(5, 2, 7, 4)}) == mm(20, 10));
}

TEST_CASE("A magnet with apex takes a single crossing shorter than its length for a touch", "[PreciseSeam][Magnet]")
{
    // A round modifier of radius 1 mm reaches 0.03 mm into the right side of the notch, square to it at
    // (22.75, 15.5): one crossing about 0.48 mm long. A vertex of the circle points straight at the side.
    const Vec2d normal(-2. / std::sqrt(5.), 1. / std::sqrt(5.)); // From the side into the notch.
    const Vec2d centre = Vec2d(22.75, 15.5) + 0.97 * normal;
    const double first = std::atan2(-normal.y(), -normal.x());
    Polygon      round;
    for (int i = 0; i < 64; ++i)
        round.points.push_back(mm(centre.x() + std::cos(first + 2 * PI * i / 64), centre.y() + std::sin(first + 2 * PI * i / 64)));
    SeamFixture                     fixture;
    PreciseSeam::EffectiveModifiers effective;
    const auto seam = [&](double length) {
        const ModelVolume *volume = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {round});
        PreciseSeam::EffectiveModifier rule = magnet_rule(ModelVolumeType::PRECISE_SEAM_CENTER);
        rule.magnet_apex        = true;
        rule.magnet_apex_length = scale_(length);
        effective.emplace(volume, rule);
        Polygon perimeter = notched();
        return PreciseSeam::insert_strong_seam_point_by_rules({volume}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer,
                                                     fixture.cache, nullptr, effective);
    };
    // Longer than 0.4 mm: the type acts, in the middle of the crossing.
    const std::optional<Point> typed = seam(0.4);
    REQUIRE(typed.has_value());
    CHECK((unscaled(*typed) - Vec2d(22.75, 15.5)).norm() < 0.01);
    // Shorter than 1 mm, or than 0.6: a touch, and the magnet finds the tip from it.
    CHECK(seam(1.) == mm(20, 10));
    CHECK(seam(0.6) == mm(20, 10));
}

TEST_CASE("A magnet with apex that touches a corner over its tip takes the tip, not an end of the touch", "[PreciseSeam][Magnet]")
{
    // A round modifier of radius 2 mm outside the corner (20, 20) of the body, its centre 1.84 mm away: the
    // corner pokes 0.16 mm into it, one crossing about 0.44 mm long. The circle is symmetric to the diagonal.
    Polygon round;
    for (int i = 0; i < 64; ++i)
        round.points.push_back(mm(21.3 + 2 * std::cos(PI / 4 + 2 * PI * i / 64), 21.3 + 2 * std::sin(PI / 4 + 2 * PI * i / 64)));
    SeamFixture                     fixture;
    PreciseSeam::EffectiveModifiers effective;
    const auto seam = [&](ModelVolumeType type) {
        const ModelVolume             *volume = fixture.add(type, {round});
        PreciseSeam::EffectiveModifier rule   = magnet_rule(type);
        rule.magnet_apex        = true;
        rule.magnet_apex_length = scale_(1.);
        effective.emplace(volume, rule);
        Polygon perimeter = rectangle(0, 0, 20, 20);
        return PreciseSeam::insert_strong_seam_point_by_rules({volume}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer,
                                                     fixture.cache, nullptr, effective);
    };
    // Center takes the middle of the touch, the corner; the ends of the crossing lie nearest to the circle.
    CHECK(seam(ModelVolumeType::PRECISE_SEAM_CENTER) == mm(20, 20));
    // Left and Right take its start on the right face and its end on the top face, as their types do; the
    // corner inside the modifier lies nearer to its centroid, so apex keeps them.
    const std::optional<Point> left  = seam(ModelVolumeType::PRECISE_SEAM_LEFT);
    const std::optional<Point> right = seam(ModelVolumeType::PRECISE_SEAM_RIGHT);
    REQUIRE(left.has_value());
    REQUIRE(right.has_value());
    CHECK((unscaled(*left) - Vec2d(20, 19.78)).norm() < 0.02);
    CHECK((unscaled(*right) - Vec2d(19.78, 20)).norm() < 0.02);
}

TEST_CASE("Strong and weak magnets with apex act on the tip between two crossings", "[PreciseSeam][Magnet]")
{
    SeamFixture                     fixture;
    PreciseSeam::EffectiveModifiers effective;
    const auto add = [&](ModelVolumeType type, bool apex) {
        const ModelVolume *volume = fixture.add(type, {rectangle(17, 13, 23, 15)});
        PreciseSeam::EffectiveModifier rule = magnet_rule(type);
        rule.magnet_apex = apex;
        effective.emplace(volume, rule);
        return volume;
    };
    const auto strong = [&](const ModelVolume *volume) {
        Polygon perimeter = notched();
        return PreciseSeam::insert_strong_seam_point_by_rules({volume}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer,
                                                     fixture.cache, nullptr, effective);
    };
    CHECK(strong(add(ModelVolumeType::PRECISE_SEAM_CENTER, true)) == mm(20, 10));
    // Without apex the crossing makes the type act: the middle of a crossing, not the tip.
    const std::optional<Point> plain = strong(add(ModelVolumeType::PRECISE_SEAM_CENTER, false));
    REQUIRE(plain.has_value());
    CHECK(*plain != mm(20, 10));
    // A weak one types the tip alone.
    const ModelVolume *green_apex = add(ModelVolumeType::PRECISE_SEAM_ENFORCED, true);
    Polygon            perimeter  = notched();
    const auto segments = PreciseSeam::collect_weak_modifier_segments_by_rules({green_apex}, perimeter, PreciseSeam::PreparedPerimeter(perimeter),
                                                                      fixture.layer, fixture.cache, nullptr, effective);
    REQUIRE(segments.size() == 1);
    CHECK(segments.front().left_point == mm(20, 10));
    CHECK(segments.front().right_point == mm(20, 10));
}

namespace {
// The one Center crossing of a perimeter with a modifier loop.
PreciseSeam::PerimeterSegment center_crossing(const Polygon &perimeter, const Polygon &loop)
{
    const PreciseSeam::PreparedPerimeter prepared(perimeter);
    const PreciseSeam::ModifierRegions   modifier  = PreciseSeam::prepare_modifier_regions(ExPolygons{ExPolygon(loop)});
    PreciseSeam::SegmentExtraction       extracted = PreciseSeam::extract_perimeter_segments(prepared, modifier,
                                                                                             ModelVolumeType::PRECISE_SEAM_CENTER);
    REQUIRE(extracted.segments.size() == 1);
    REQUIRE(extracted.segments.front().strong_target.has_value());
    return std::move(extracted.segments.front());
}

// Where apex moves the middle of that crossing, with a threshold in mm.
Point segment_apex(const Polygon &perimeter, const Polygon &loop, double threshold)
{
    const PreciseSeam::PerimeterSegment segment = center_crossing(perimeter, loop);
    return PreciseSeam::find_segment_apex(segment, *segment.strong_target, scale_(threshold)).point;
}
} // namespace

TEST_CASE("apex moves the middle of a Center crossing to the tip farthest from its chord", "[PreciseSeam][Apex]")
{
    const size_t shift = GENERATE(0, 3, 5);
    CAPTURE(shift);
    SECTION("A tip deeper than the middle by more than the threshold")
    {
        // The crossing runs (23, 16) - (20, 10) - (19, 12); its middle (21, 12) lies 1.414 mm from the chord, the
        // tip 2.121 mm: 0.707 mm deeper.
        const Polygon perimeter = rotated(notched(), shift);
        const Polygon loop      = rectangle(19, 8, 30, 16);
        CHECK((unscaled(center_crossing(perimeter, loop).strong_target->point) - Vec2d(21, 12)).norm() < 0.001);
        CHECK(segment_apex(perimeter, loop, 0.001) == mm(20, 10));
        CHECK(segment_apex(perimeter, loop, 0.7) == mm(20, 10));
        CHECK(segment_apex(perimeter, loop, 0.71) != mm(20, 10));
    }
    SECTION("A middle on the chord, or a crossing without inner vertices, keeps its place")
    {
        // Over a corner symmetrically: the middle is the corner itself, as deep as any vertex.
        const Polygon square = rotated(rectangle(0, 0, 20, 20), shift % 4);
        CHECK(segment_apex(square, rectangle(18, 18, 22, 22), 0.001) == mm(20, 20));
        // Across one face: two points only.
        const PreciseSeam::PerimeterSegment face = center_crossing(square, rectangle(8, -2, 12, 2));
        CHECK(PreciseSeam::find_segment_apex(face, *face.strong_target, scale_(0.001)).point == face.strong_target->point);
    }
    SECTION("A vertex on the other side of the chord is skipped, the others still count")
    {
        // A step in the front face: the crossing runs (20, 0) - (26, 0) - (26, 2) - (28, 2), its middle (25, 0).
        // (26, 0) lies on the side of the middle and deeper; (26, 2) on the other side.
        const Polygon stepped(Points{mm(0, 0), mm(26, 0), mm(26, 2), mm(40, 2), mm(40, 20), mm(0, 20)});
        CHECK(segment_apex(rotated(stepped, shift), rectangle(20, -1, 28, 3), 0.001) == mm(26, 0));
    }
    SECTION("Equally deep tips go to the one nearest the middle")
    {
        // Two V notches 8 mm deep under a horizontal chord; the ridge between them is off centre, and the middle
        // (22.65, 15.58) lies nearer the right tip, which rear-then-left alone would not pick.
        const Polygon twin(Points{mm(0, 0), mm(40, 0), mm(40, 20), mm(28, 20), mm(26, 10), mm(23, 16), mm(18, 10), mm(16, 20), mm(0, 20)});
        CHECK(segment_apex(rotated(twin, shift), rectangle(15, 5, 30, 18), 0.001) == mm(26, 10));
    }
}

TEST_CASE("The strong driver applies apex to Center only, with the print resolution by default", "[PreciseSeam][Apex]")
{
    SeamFixture                     fixture;
    PreciseSeam::EffectiveModifiers effective;
    const auto seam = [&](ModelVolumeType type, bool apex, std::optional<double> threshold = std::nullopt) {
        const ModelVolume             *volume = fixture.add(type, {rectangle(19, 8, 30, 16)});
        PreciseSeam::EffectiveModifier rule{type, {}};
        rule.apex = apex;
        if (threshold)
            rule.apex_threshold = scale_(*threshold);
        effective.emplace(volume, rule);
        Polygon perimeter = notched();
        return PreciseSeam::insert_strong_seam_point_by_rules({volume}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer,
                                                              fixture.cache, nullptr, effective);
    };
    REQUIRE(fixture.print.config().resolution.value < 0.7);
    CHECK(seam(ModelVolumeType::PRECISE_SEAM_CENTER, true) == mm(20, 10));
    CHECK(seam(ModelVolumeType::PRECISE_SEAM_CENTER, true, 0.5) == mm(20, 10));
    const std::optional<Point> middle = seam(ModelVolumeType::PRECISE_SEAM_CENTER, false);
    REQUIRE(middle.has_value());
    CHECK((unscaled(*middle) - Vec2d(21, 12)).norm() < 0.001);
    CHECK(seam(ModelVolumeType::PRECISE_SEAM_CENTER, true, 1.) == middle);
    // Left keeps the start of the crossing.
    const std::optional<Point> left = seam(ModelVolumeType::PRECISE_SEAM_LEFT, true);
    REQUIRE(left.has_value());
    CHECK((unscaled(*left) - Vec2d(23, 16)).norm() < 0.001);
}

TEST_CASE("A strong magnet places the seam where it does not cross the perimeter, and acts by its type where it does", "[PreciseSeam][Magnet]")
{
    SeamFixture                     fixture;
    PreciseSeam::PreciseSeamWarnings warnings;
    PreciseSeam::EffectiveModifiers effective;
    const auto add = [&](ModelVolumeType type, Polygon slice, std::optional<double> limit = std::nullopt) {
        const ModelVolume *volume = fixture.add(type, {std::move(slice)});
        effective.emplace(volume, magnet_rule(type, limit));
        warnings.modifier_usage.try_emplace(volume);
        return volume;
    };
    const auto seam = [&](const std::vector<const ModelVolume*> &volumes) {
        Polygon perimeter = rectangle(0, 0, 20, 20);
        return PreciseSeam::insert_strong_seam_point_by_rules(volumes, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer,
                                                     fixture.cache, &warnings, effective);
    };
    // Left is the start of the crossing on the counterclockwise front face, and means nothing to a magnet.
    CHECK(seam({add(ModelVolumeType::PRECISE_SEAM_LEFT, rectangle(8, -5, 12, -3))}) == mm(10, 0));
    CHECK(seam({add(ModelVolumeType::PRECISE_SEAM_LEFT, rectangle(8, -2, 12, 2))}) == mm(8, 0));
    // Full containment is no crossing: the nearest point instead of a skip, without its warning.
    const ModelVolume *around = add(ModelVolumeType::PRECISE_SEAM_CENTER, rectangle(-3, -2, 25, 30));
    CHECK(seam({around}) == mm(11, 0));
    CHECK(warnings.full_containment.load() == 0);
    CHECK(warnings.modifier_usage.at(around).reached.load());
    // Beyond its limit the magnet passes the perimeter on, and has not reached it.
    const ModelVolume *far = add(ModelVolumeType::PRECISE_SEAM_CENTER, rectangle(8, -5, 12, -3), 2.);
    const ModelVolume *junior = fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, {rectangle(14, -2, 16, 2)});
    CHECK(seam({far, junior}) == mm(15, 0));
    CHECK(warnings.modifier_usage.at(far).checked.load());
    CHECK_FALSE(warnings.modifier_usage.at(far).reached.load());
    // A limited ring around the body is too far as well: full containment does not count as reached.
    const ModelVolume *wide = add(ModelVolumeType::PRECISE_SEAM_CENTER, rectangle(-5, -5, 25, 25), 2.);
    CHECK_FALSE(seam({wide}).has_value());
    CHECK_FALSE(warnings.modifier_usage.at(wide).reached.load());
    // A senior magnet takes the perimeter from a junior ordinary modifier.
    CHECK(seam({add(ModelVolumeType::PRECISE_SEAM_CENTER, rectangle(8, -5, 12, -3)), junior}) == mm(10, 0));
}

TEST_CASE("A weak magnet types its point only", "[PreciseSeam][Magnet]")
{
    using Seam = SeamPlacerImpl::EnforcedBlockedSeamPoint;
    const auto [type, painted, typed] = GENERATE(table<ModelVolumeType, Seam, Seam>({
        {ModelVolumeType::PRECISE_SEAM_ENFORCED, Seam::Blocked, Seam::Enforced},
        {ModelVolumeType::PRECISE_SEAM_BLOCKED, Seam::Enforced, Seam::Blocked},
        {ModelVolumeType::PRECISE_SEAM_NEUTRAL, Seam::Enforced, Seam::Neutral}}));
    // In front of the body, and around it: full containment is no crossing, so one point either way.
    const auto [slice, point] = GENERATE(table<Polygon, Point>({{rectangle(8, -5, 12, -3), mm(10, 0)}, {rectangle(-3, -2, 25, 30), mm(11, 0)}}));
    SeamFixture                     fixture;
    PreciseSeam::EffectiveModifiers effective;
    const ModelVolume              *volume = fixture.add(type, {slice});
    effective.emplace(volume, magnet_rule(type));
    Polygon    perimeter = rectangle(0, 0, 20, 20);
    const auto segments  = PreciseSeam::collect_weak_modifier_segments_by_rules({volume}, perimeter, PreciseSeam::PreparedPerimeter(perimeter),
                                                                       fixture.layer, fixture.cache, nullptr, effective);
    REQUIRE(segments.size() == 1);
    CHECK_FALSE(segments.front().whole_perimeter);
    const auto types = weak_candidate_types(perimeter, segments, [painted = painted](const Point &) { return painted; });
    REQUIRE(types.size() == perimeter.size());
    for (size_t i = 0; i < perimeter.size(); ++i) {
        CAPTURE(i, perimeter[i].x(), perimeter[i].y());
        CHECK(types[i] == (perimeter[i] == point ? typed : painted));
    }
}

TEST_CASE("The rules give a magnet its limit in scaled units, a touch and no limit apart", "[PreciseSeam][Rules][Magnet]")
{
    ExportFixture      fixture("back");
    const ModelVolume *any   = add_modifier(fixture, green, "magnet");
    const ModelVolume *touch = add_modifier(fixture, green, "magnet0");
    const ModelVolume *tiny  = add_modifier(fixture, green, "magnet0.000001");
    const ModelVolume *plain = add_modifier(fixture, green, "");
    const ModelVolume *apex  = add_modifier(fixture, green, "magnet2apex0.5");
    const WeakInit     init(fixture, 10);
    CHECK(init.effective.at(any).magnet);
    CHECK_FALSE(init.effective.at(any).magnet_limit.has_value());
    CHECK_FALSE(init.effective.at(any).magnet_apex);
    CHECK(init.effective.at(apex).magnet_apex);
    CHECK(init.effective.at(apex).magnet_limit == scale_(2.));
    CHECK(init.effective.at(apex).magnet_apex_length == scale_(0.5));
    CHECK(init.effective.at(any).magnet_apex_length == 0.);
    CHECK(init.effective.at(touch).magnet_limit == 0.);
    // One nanometre stays a positive limit, not a touch and not unlimited: one coordinate unit at the default
    // scale, a tenth of one on large beds, where a whole-unit limit would round it to zero.
    REQUIRE(init.effective.at(tiny).magnet_limit.has_value());
    CHECK_THAT(*init.effective.at(tiny).magnet_limit, Catch::Matchers::WithinRel(scale_(0.000001), 1e-9));
    CHECK(*init.effective.at(tiny).magnet_limit > 0.);
    CHECK_FALSE(init.effective.at(plain).magnet);
}

TEST_CASE("A magnet modifier draws the exported seam to the nearest point of the wall", "[PreciseSeam][Rules][Magnet]")
{
    // The helper stands 3 mm in front of the front face, around x = 10 mm: the seam goes to the face centre.
    const auto [type, rules, drawn] = GENERATE(table<ModelVolumeType, std::string, bool>({
        {ModelVolumeType::PRECISE_SEAM_CENTER, "magnet", true},
        {ModelVolumeType::PRECISE_SEAM_CENTER, "magnet5", true},
        {ModelVolumeType::PRECISE_SEAM_CENTER, "magnet2", false},
        {ModelVolumeType::PRECISE_SEAM_ENFORCED, "magnet", true}}));
    CAPTURE(rules);
    ExportFixture fixture("back");
    fixture.add(type, make_cube(2, 2, 4), Vec3d(9, -5, -1));
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = rules;
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (const OuterLoop &loop : loops) {
        CAPTURE(loop.seam.x(), loop.seam.y());
        if (drawn)
            CHECK(front_helper(loop, {0.}) == 0);
        else
            CHECK(seam_at_back(loop));
    }
    // A magnet beyond its limit on every layer had no effect, which is reported.
    CHECK(fixture.warning().empty() == drawn);
}

TEST_CASE("apex puts the exported seam of a Center modifier into the corner it crosses", "[PreciseSeam][Rules][Apex]")
{
    // The helper covers the front left corner unevenly: the outer wall crosses it 0.8 mm down the left face,
    // then 2.8 mm along the front face, so the middle lies 1 mm right of the corner (0.2, 0.2), and the corner
    // 0.275 mm farther from the chord than the middle.
    const auto [type, rules, x, y, ignored] = GENERATE(table<ModelVolumeType, std::string, double, double, bool>({
        {ModelVolumeType::PRECISE_SEAM_CENTER, "", 1.2, 0.2, false},
        {ModelVolumeType::PRECISE_SEAM_CENTER, "apex", 0.2, 0.2, false},
        {ModelVolumeType::PRECISE_SEAM_LEFT, "center apex", 0.2, 0.2, false},
        {ModelVolumeType::PRECISE_SEAM_CENTER, "apex0.25", 0.2, 0.2, false},
        {ModelVolumeType::PRECISE_SEAM_CENTER, "apex0.3", 1.2, 0.2, false},
        // Left takes the start of the crossing, on the left face; its apex is reported and ignored.
        {ModelVolumeType::PRECISE_SEAM_LEFT, "apex", 0.2, 1., true},
        {ModelVolumeType::PRECISE_SEAM_CENTER, "left apex", 0.2, 1., true}}));
    CAPTURE(rules);
    ExportFixture fixture("back");
    fixture.add(type, make_cube(4, 2, 4), Vec3d(-1, -1, -1));
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = rules;
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    for (const OuterLoop &loop : loops) {
        CAPTURE(loop.seam.x(), loop.seam.y());
        CHECK_THAT(loop.seam.x(), Catch::Matchers::WithinAbs(x - 10., 0.01));
        CHECK_THAT(loop.seam.y(), Catch::Matchers::WithinAbs(y - 10., 0.01));
    }
    CHECK((fixture.warning().find("apex acts only with the Seam Center type") != std::string::npos) == ignored);
}

TEST_CASE("Magnet relay members take the seam in turn", "[PreciseSeam][Rules][Relay][Magnet]")
{
    // A and B stand 3 mm in front around x = 5 and 15 mm; B's limit of 2 mm keeps it from every layer.
    const bool limited = GENERATE(false, true);
    CAPTURE(limited);
    ExportFixture fixture("back");
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 2, 4), Vec3d(4, -5, -1), "A");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = "magnet,relay";
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 2, 4), Vec3d(14, -5, -1), "B");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = limited ? "magnet2,relay" : "magnet,relay";
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    // Where B cannot, A stands in on every layer.
    CHECK(front_helpers(loops, {-5., 5.}) ==
          (limited ? std::vector<int>(10, 0) : std::vector<int>{0, 1, 0, 1, 0, 1, 0, 1, 0, 1}));
}

// Expected holders come from _myDocs/_relay/26.10.24_ps-seed-ui/relay_seed_oracle.py.
TEST_CASE("A strong relay-random group takes the seeds of its members, one without layers included", "[PreciseSeam][Rules][Relay]")
{
    // A, B and C as in the relay-random test above; D follows C in the queue, and its rules exclude every layer.
    const auto [seeds, expected] = GENERATE(table<std::vector<std::string>, std::vector<int>>({
        {{"", "", "", ""}, {0, 1, 2, 1, 2, 0, 2, 1, 2, 0}},
        {{"seed 7", "", "", ""}, {1, 0, 2, 0, 1, 0, 1, 2, 0, 2}},
        {{"", "", "", "seed 7"}, {0, 1, 2, 1, 0, 1, 0, 1, 2, 1}}}));
    CAPTURE(seeds);
    ExportFixture fixture("back");
    const std::vector<double> xs{2., 9., 16., 9.};
    for (size_t i = 0; i < xs.size(); ++i) {
        // B is Left, so its seam lies on its left edge, as in the test above.
        fixture.add(i == 1 ? ModelVolumeType::PRECISE_SEAM_LEFT : ModelVolumeType::PRECISE_SEAM_CENTER, front_box(), Vec3d(xs[i], -3, -1),
                    std::string(1, char('A' + i)));
        fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random " + seeds[i] + (i == 3 ? ", 500" : "");
    }
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_helpers(loops, {-7., -1., 7.}) == expected);
}

TEST_CASE("A weak relay-random group takes the seeds of its members", "[PreciseSeam][Rules][Relay]")
{
    const auto [seeds, expected] = GENERATE(table<std::vector<std::string>, std::vector<int>>({
        {{"", "", "", ""}, {0, 1, 2, 1, 2, 0, 2, 1, 2, 0}},
        {{"", "seed 7", "", ""}, {0, 2, 0, 1, 0, 2, 1, 2, 0, 1}},
        {{"", "", "", "seed 7"}, {0, 1, 2, 1, 0, 1, 0, 1, 2, 1}}}));
    CAPTURE(seeds);
    ExportFixture fixture("back");
    const std::vector<double> xs{2., 9., 16., 9.};
    for (size_t i = 0; i < xs.size(); ++i) {
        fixture.add(ModelVolumeType::PRECISE_SEAM_ENFORCED, front_box(), Vec3d(xs[i], -3, -1));
        fixture.model.objects.front()->volumes.back()->precise_seam_rules = "relay-random " + seeds[i] + (i == 3 ? ", 500" : "");
    }
    const auto loops = outer_wall_loops(fixture.export_gcode());
    REQUIRE(loops.size() == 10);
    CHECK(front_zones(loops, {2., 9., 16.}) == expected);
}

TEST_CASE("model counts layers and heights from the model, so a raft does not move them", "[PreciseSeam][Rules][Origin]")
{
    const int raft_layers = GENERATE(0, 2);
    const std::string rule = GENERATE(as<std::string>{}, "model, 3-5", "model, excl, 1-8", "model, z0-1", "model, zin0.4-1");
    CAPTURE(raft_layers, rule);
    ExportFixture fixture("back", raft_layers);
    fixture.add(ModelVolumeType::PRECISE_SEAM_CENTER, make_cube(2, 6, 4), Vec3d(4, -3, -1), "front");
    fixture.model.objects.front()->volumes.back()->precise_seam_rules = rule;
    const auto loops = outer_wall_loops(fixture.export_gcode());
    const PrintObject &object = *fixture.print.objects().front();
    REQUIRE(loops.size() == 10);
    REQUIRE(object.layers().size() == loops.size());
    // Loop i is model layer i + 1; its band is measured from the bottom of the first model layer.
    const double bottom = object.layers().front()->bottom_z();
    for (size_t i = 0; i < loops.size(); ++i) {
        const double low = object.layers()[i]->bottom_z() - bottom, high = object.layers()[i]->print_z - bottom;
        CAPTURE(i, low, high);
        const bool expected = rule == "model, 3-5"       ? i >= 2 && i <= 4 :
                              rule == "model, excl, 1-8" ? i >= 8 :
                              rule == "model, z0-1"      ? low < 1. - EPSILON :
                                                           low >= 0.4 - EPSILON && high <= 1. + EPSILON;
        CHECK(seam_on_front_helper(loops[i]) == expected);
    }
    CHECK(fixture.warning().empty());
}
