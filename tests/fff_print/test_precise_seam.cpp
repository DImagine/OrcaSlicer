#include <catch2/catch_all.hpp>

#include "test_helpers.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/GCode/PreciseSeam.hpp"

#include <algorithm>

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
    const Polygon &polygon, const std::vector<PreciseSeam::WeakModifierSegment> &segments)
{
    // Use the same coordinate conversion as production when applying prepared boundaries.
    PrintObjectSeamData::LayerSeams candidates;
    candidates.perimeters.emplace_back();
    auto &loop = candidates.perimeters.back();
    loop.start_index = 0;
    for (const Point &point : polygon.points) {
        const Vec2f position = unscale(point).cast<float>();
        candidates.points.emplace_back(Vec3f(position.x(), position.y(), 0), loop, 0,
                                       SeamPlacerImpl::EnforcedBlockedSeamPoint::Neutral);
    }
    loop.end_index = candidates.points.size();
    bool enforced = false;
    PreciseSeam::apply_weak_modifiers_to_perimeter(segments, candidates, loop, enforced);
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
    // Compare with the shared flattened slicing path used by support modifiers.
    const auto flat_layers = object->slice_modifier_volumes({volume});
    REQUIRE(flat_layers.size() == 1);
    CHECK(flat_layers.front().size() == 3);
    double flat_area = 0.;
    for (const Polygon &contour : flat_layers.front())
        flat_area += contour.area();
    CHECK_THAT(flat_area, Catch::Matchers::WithinAbs(area, 1.));
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
    CHECK(warnings.multiple_intersections.load() == extra_segment);
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
    // Duplicate boundaries used to stall the reverse cursor before reaching the separate segment.
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

TEST_CASE("Weak full containment keeps its warning and leaves the perimeter unchanged", "[PreciseSeam]")
{
    const auto type = GENERATE(ModelVolumeType::PRECISE_SEAM_ENFORCED, ModelVolumeType::PRECISE_SEAM_BLOCKED,
                               ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    SeamFixture fixture;
    Polygon perimeter = rectangle(0, 0, 20, 20);
    const Points original = perimeter.points;
    const auto *modifier = fixture.add(type, {rectangle(-2, -2, 22, 22)});
    PreciseSeam::PreciseSeamWarnings warnings;
    const auto segments = PreciseSeam::collect_weak_modifier_segments(
        {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings);
    CHECK(segments.empty());
    CHECK(perimeter.points == original);
    CHECK(warnings.full_containment.load());
    CHECK_FALSE(warnings.multiple_intersections.load());
    CHECK(warnings.failed_fragments.load() == 0);
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
    else
        CHECK(PreciseSeam::collect_weak_modifier_segments(
            {modifier}, perimeter, PreciseSeam::PreparedPerimeter(perimeter), fixture.layer, fixture.cache, &warnings).empty());
    CHECK(perimeter.points == original);
    CHECK(warnings.full_containment.load());
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
    CHECK(warnings.full_containment.load() == (scenario == 2 || scenario == 3));
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
    CHECK(warnings.full_containment.load() == contained);
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
    CHECK(warnings.multiple_intersections.load() == strong); // Warn only strong, after collecting all segments.
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
