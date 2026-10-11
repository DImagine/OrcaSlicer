#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/GCode/SeamChain.hpp"
#include "libslic3r/GCode/SeamFeature.hpp"
#include "libslic3r/Point.hpp"
#include <admesh/stl.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::SeamPlacerImpl;

namespace {
std::vector<SeamChainPoint> straight_chain(size_t count = 21, double size = 5.)
{
    std::vector<SeamChainPoint> chain;
    for (size_t i = 0; i < count; ++i) {
        SeamChainPoint p;
        p.width = .4;
        p.vertices = {{Vec2d(0., 0.)}, {Vec2d(size, 0.)}, {Vec2d(size, 5.)}, {Vec2d(0., 5.)}};
        p.initial = Vec3d(.1, 0., .3 + .1 * double(i));
        p.candidate_arc = p.anchor_arc = .1;
        chain.push_back(p);
    }
    return chain;
}

indexed_triangle_set box_mesh(double size = 5.)
{
    indexed_triangle_set mesh;
    for (double z : {0., 3.}) for (double y : {-.2, 5.2}) for (double x : {-.2, size + .2})
        mesh.vertices.emplace_back(float(x), float(y), float(z));
    const int quads[6][4] = {{0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4},
                             {2, 3, 7, 6}, {0, 2, 6, 4}, {1, 3, 7, 5}};
    const Vec3f center(float(size / 2.), 2.5f, 1.5f);
    for (const auto &q : quads) for (int j : {1, 2}) {
        Vec3i32 face(q[0], q[j], q[j + 1]);
        const Vec3f a = mesh.vertices[size_t(face[0])], b = mesh.vertices[size_t(face[1])], c = mesh.vertices[size_t(face[2])];
        if ((b - a).cross(c - a).dot((a + b + c) / 3.f - center) < 0.) std::swap(face[1], face[2]);
        mesh.indices.push_back(face); // Fixtures have explicit outward winding and shared vertex identities.
    }
    return mesh;
}

std::vector<SeamFeatureBinding> middle_targets(const std::vector<SeamChainPoint> &chain)
{
    std::vector<SeamFeatureBinding> targets(chain.size());
    for (size_t i = 8; i <= 12; ++i)
        targets[i] = {SeamFeatureState::Proposed, 0, 0, Vec3d(0., 0., chain[i].initial.z()), 0.};
    return targets;
}
} // namespace

TEST_CASE("Feature targets remain exact while free neighbors form a bounded transition", "[SeamPlacer][SeamFeature]")
{
    auto chain = straight_chain();
    chain.front().fixed = true; // Represents an individual strong Precise Seam, not an excluded whole chain.
    const auto result = accept_seam_feature_targets(chain, middle_targets(chain));
    REQUIRE(result.accepted);
    for (size_t i = 8; i <= 12; ++i) {
        CHECK(result.bindings[i].state == SeamFeatureState::Locked);
        CHECK_THAT(result.positions[i].x(), Catch::Matchers::WithinAbs(0., 1e-12));
        CHECK(result.locked[i]);
    }
    CHECK((result.positions.front() - chain.front().initial).norm() < 1e-12);
    CHECK(result.positions[7].x() < chain[7].initial.x());
    CHECK_FALSE(result.locked[7]);
    CHECK(seam_feature_continuous(chain, result.positions, result.geometric));
    for (size_t i = 0; i < chain.size(); ++i) {
        chain[i].initial = result.positions[i];
        chain[i].fixed = result.locked[i];
    }
    SeamChainOptions options;
    options.independent_intervals = true;
    const auto optimized = optimize_seam_chain(chain, options);
    for (size_t i = 0; i < chain.size(); ++i) if (result.locked[i])
        CHECK((optimized.positions[i] - result.positions[i]).norm() < 1e-12);
}

TEST_CASE("A feature target outside its painted interval is rejected before freezing", "[SeamPlacer][SeamFeature]")
{
    auto chain = straight_chain();
    auto targets = middle_targets(chain);
    for (size_t i = 8; i <= 12; ++i) chain[i].lower_arc = .05;
    const auto result = accept_seam_feature_targets(chain, targets);
    CHECK_FALSE(result.accepted);
    for (size_t i = 0; i < chain.size(); ++i) CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
}

TEST_CASE("Infeasible transition packets retain every original point", "[SeamPlacer][SeamFeature][Regression]")
{
    auto chain = straight_chain();
    for (size_t i = 0; i < chain.size(); ++i) chain[i].fixed = i < 8 || i > 12;
    SeamFeatureOptions options;
    options.step_growth = options.step_allowance_widths = options.turn_allowance = 0.;
    const auto result = accept_seam_feature_targets(chain, middle_targets(chain), options);
    CHECK_FALSE(result.accepted);
    for (size_t i = 0; i < chain.size(); ++i) {
        CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
        CHECK(result.locked[i] == chain[i].fixed);
    }
    CHECK(result.bindings[10].state == SeamFeatureState::TransitionRejected);
}

TEST_CASE("An exhausted transition solver commits no provisional anchors", "[SeamPlacer][SeamFeature]")
{
    auto chain = straight_chain();
    for (size_t i = 0; i < chain.size(); ++i) if (i < 8 || i > 12) chain[i].lower_arc = .095;
    SeamFeatureOptions options;
    options.max_iterations = 1;
    const auto result = accept_seam_feature_targets(chain, middle_targets(chain), options);
    CHECK_FALSE(result.accepted);
    for (size_t i = 0; i < chain.size(); ++i) CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
}

TEST_CASE("Mesh feature discovery binds the inset printed corner rather than the mesh edge", "[SeamPlacer][SeamFeature]")
{
    const auto chain = straight_chain();
    const auto results = bind_seam_features(box_mesh(), {chain});
    REQUIRE(results.size() == 1);
    REQUIRE(results[0].accepted);
    for (size_t i = 0; i < chain.size(); ++i) {
        CHECK_THAT(results[0].positions[i].x(), Catch::Matchers::WithinAbs(0., 1e-12));
        CHECK_THAT(results[0].positions[i].y(), Catch::Matchers::WithinAbs(0., 1e-12));
        CHECK(results[0].locked[i]);
    }
}

TEST_CASE("Two plausible nearby creases remain ambiguous instead of choosing the nearer one", "[SeamPlacer][SeamFeature]")
{
    auto chain = straight_chain(21, .3);
    for (auto &p : chain) { p.initial.x() = .14; p.candidate_arc = p.anchor_arc = .14; }
    const auto results = bind_seam_features(box_mesh(.3), {chain});
    REQUIRE(results.size() == 1);
    CHECK_FALSE(results[0].accepted);
    for (size_t i = 0; i < chain.size(); ++i) {
        CHECK(results[0].bindings[i].state == SeamFeatureState::Ambiguous);
        CHECK((results[0].positions[i] - chain[i].initial).norm() < 1e-12);
    }
}

TEST_CASE("Feature discovery does not depend on the order of overlapping chains", "[SeamPlacer][SeamFeature]")
{
    auto a = straight_chain(), b = straight_chain();
    for (auto &p : b) { p.initial.x() = .12; p.candidate_arc = p.anchor_arc = .12; }
    const auto first = bind_seam_features(box_mesh(), {a, b});
    const auto reverse = bind_seam_features(box_mesh(), {b, a});
    REQUIRE(first[0].accepted);
    REQUIRE(reverse[1].accepted);
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK(first[0].bindings[i].edge == reverse[1].bindings[i].edge);
        CHECK((first[0].positions[i] - reverse[1].positions[i]).norm() < 1e-12);
    }
}

TEST_CASE("Insufficient layer evidence and a search budget cannot create mesh locks", "[SeamPlacer][SeamFeature]")
{
    const auto short_chain = straight_chain(2);
    const auto short_result = bind_seam_features(box_mesh(), {short_chain});
    CHECK_FALSE(short_result[0].accepted);
    CHECK(short_result[0].bindings[0].state == SeamFeatureState::InsufficientEvidence);
    SeamFeatureOptions options;
    options.max_edge_queries = 1;
    const auto result = bind_seam_features(box_mesh(), {straight_chain()}, options);
    CHECK_FALSE(result[0].accepted);
    CHECK(result[0].bindings[0].state == SeamFeatureState::Budget);
}

TEST_CASE("Blocked contour corners and unsupported targets never become feature locks", "[SeamPlacer][SeamFeature]")
{
    auto chain = straight_chain();
    for (auto &p : chain) p.vertices[0].feature_allowed = false;
    CHECK_FALSE(bind_seam_features(box_mesh(), {chain})[0].accepted);
    chain = straight_chain();
    // The old point is 0.1 mm from this corner; make its support loss exceed the 0.1-width allowance.
    for (auto &p : chain) p.vertices[0].overhang = 10.;
    CHECK_FALSE(bind_seam_features(box_mesh(), {chain})[0].accepted);
}

TEST_CASE("Confirmed feature junctions keep their bend while artificial transition kinks are rejected", "[SeamPlacer][SeamFeature]")
{
    const auto chain = straight_chain(3);
    std::vector<Vec3d> positions;
    for (const auto &p : chain) positions.push_back(p.initial);
    positions[1].x() += .1; // Equal-size steps can still hide a visible new kink.
    CHECK_FALSE(seam_feature_continuous(chain, positions, std::vector<bool>(3, false)));
    CHECK(seam_feature_continuous(chain, positions, std::vector<bool>(3, true)));
}

TEST_CASE("Independent free intervals retain degrees of freedom beside many hard rows", "[SeamPlacer][SeamChain][SeamFeature]")
{
    auto chain = straight_chain(41);
    for (size_t i = 0; i < chain.size(); ++i) {
        chain[i].initial.x() += .03 * std::sin(double(i) * .5);
        chain[i].candidate_arc = chain[i].anchor_arc = chain[i].initial.x();
        chain[i].fixed = i <= 20;
    }
    SeamChainOptions options;
    options.segments = 1;
    CHECK(optimize_seam_chain(chain, options).status == SeamChainStatus::NoFreedom);
    options.independent_intervals = true;
    const auto result = optimize_seam_chain(chain, options);
    CHECK(result.status != SeamChainStatus::NoFreedom);
    for (size_t i = 0; i <= 20; ++i) CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
}

TEST_CASE("Incident carriers at a four-edge junction certify one printed target", "[SeamPlacer][SeamFeature]")
{
    indexed_triangle_set mesh;
    mesh.vertices = {Vec3f(-3.f, 0.f, 3.f), Vec3f(0.f, -3.f, 3.f), Vec3f(3.f, 0.f, 3.f),
                     Vec3f(0.f, 3.f, 3.f), Vec3f(0.f, 0.f, 0.f), Vec3f(0.f, 0.f, 6.f)};
    for (int apex : {4, 5}) for (int i = 0; i < 4; ++i) {
        Vec3i32 face(i, (i + 1) % 4, apex);
        const Vec3f a = mesh.vertices[size_t(face[0])], b = mesh.vertices[size_t(face[1])], c = mesh.vertices[size_t(face[2])];
        if ((b - a).cross(c - a).dot((a + b + c) / 3.f - Vec3f(0.f, 0.f, 3.f)) < 0.f) std::swap(face[1], face[2]);
        mesh.indices.push_back(face);
    }
    std::vector<SeamChainPoint> chain;
    for (size_t i = 0; i < 23; ++i) {
        const double z = .8 + .2 * double(i);
        const double radius = std::min(z, 6. - z) - .2 * std::sqrt(2.);
        SeamChainPoint p;
        p.width = .4;
        p.vertices = {{Vec2d(-radius, 0.)}, {Vec2d(0., -radius)}, {Vec2d(radius, 0.)}, {Vec2d(0., radius)}};
        p.initial = Vec3d(-radius + .1, -.1, z);
        p.candidate_arc = p.anchor_arc = .1 * std::sqrt(2.);
        chain.push_back(p);
    }
    const auto result = bind_seam_features(mesh, {chain})[0];
    REQUIRE(result.accepted);
    const size_t junction = 11;
    CHECK(result.bindings[junction].state == SeamFeatureState::Locked);
    CHECK(result.locked[junction]);
    for (size_t i = 0; i < chain.size(); ++i) if (result.locked[i])
        CHECK((result.positions[i].head<2>() - chain[i].vertices[0].position).norm() < 1e-12);
}

TEST_CASE("Coincident but disconnected mesh components never form a shared feature junction", "[SeamPlacer][SeamFeature]")
{
    auto mesh = box_mesh();
    const auto copy = mesh;
    const int offset = int(mesh.vertices.size());
    mesh.vertices.insert(mesh.vertices.end(), copy.vertices.begin(), copy.vertices.end());
    for (auto face : copy.indices) {
        face.array() += offset; // Identical coordinates deliberately retain separate topological identities.
        mesh.indices.push_back(face);
    }
    const auto result = bind_seam_features(mesh, {straight_chain()})[0];
    CHECK_FALSE(result.accepted);
    CHECK(result.bindings[10].state == SeamFeatureState::Ambiguous);
}

TEST_CASE("Invalid feature controls retain the complete initial line", "[SeamPlacer][SeamFeature]")
{
    const auto chain = straight_chain();
    SeamFeatureOptions options;
    options.arc_chart_fraction = .5; // The periodic arc chart must remain smaller than half the loop.
    const auto result = bind_seam_features(box_mesh(), {chain}, options)[0];
    CHECK_FALSE(result.accepted);
    for (size_t i = 0; i < chain.size(); ++i) CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
}

TEST_CASE("A bounded local repair finds a continuous bridge that the quadratic seed misses", "[SeamPlacer][SeamFeature][Regression]")
{
    const double z[] = {.1, .28, .44, .6, .76, .92};
    const double x[] = {.739269, .549805, .404537, .259793, .094770, -.089545};
    const double y[] = {0., -.22, -.34, -.34, -.34, -.34};
    std::vector<SeamChainPoint> points;
    std::vector<SeamFeatureBinding> targets(6);
    for (size_t i = 0; i < 6; ++i) {
        const double left = i < 3 ? 0. : -.16 * double(i - 3);
        SeamChainPoint p;
        p.width = .42;
        p.anchor_weight = 3.6976744;
        p.initial = Vec3d(x[i], y[i], z[i]);
        p.vertices = {{Vec2d(left, y[i])}, {Vec2d(left + 5., y[i])}, {Vec2d(left + 5., y[i] + 5.)}, {Vec2d(left, y[i] + 5.)}};
        p.candidate_arc = p.anchor_arc = x[i] - left;
        if (i == 1) p.lower_arc = p.candidate_arc - .09072;
        if (i == 2) p.lower_arc = p.candidate_arc - .22365;
        if (i >= 3) targets[i] = {SeamFeatureState::Proposed, 0, 0, Vec3d(left, y[i], z[i]), 0.};
        points.push_back(p);
    }
    SeamFeatureOptions options;
    options.bridge_repair_sweeps = 0;
    CHECK_FALSE(accept_seam_feature_targets(points, targets, options).accepted);
    options.bridge_repair_sweeps = 8;
    const auto result = accept_seam_feature_targets(points, targets, options);
    REQUIRE(result.accepted);
    CHECK(seam_feature_continuous(points, result.positions, result.geometric, options));
    for (size_t i = 3; i < points.size(); ++i) {
        CHECK(result.locked[i]);
        CHECK((result.positions[i] - targets[i].target).norm() < 1e-12);
    }
}

TEST_CASE("Verified points on one inclined carrier use its geometric layer spacing", "[SeamPlacer][SeamFeature][Regression]")
{
    auto chain = straight_chain(3);
    std::vector<SeamFeatureBinding> targets(3);
    const Vec3d direction = Vec3d(.08, 0., .1).normalized();
    for (size_t i = 0; i < chain.size(); ++i) {
        const double left = .08 * double(i);
        chain[i].initial.x() = .25; // The initial spline compresses lateral motion along the inclined crease.
        chain[i].vertices = {{Vec2d(left, 0.)}, {Vec2d(left + 5., 0.)}, {Vec2d(left + 5., 5.)}, {Vec2d(left, 5.)}};
        chain[i].candidate_arc = chain[i].anchor_arc = .25 - left;
        targets[i] = {SeamFeatureState::Proposed, 7, 0, Vec3d(left, 0., chain[i].initial.z()), 0., direction};
    }
    SeamFeatureOptions options;
    options.step_growth = options.step_allowance_widths = 0.;
    options.carrier_reference_steps = false;
    CHECK_FALSE(accept_seam_feature_targets(chain, targets, options).accepted);
    options.carrier_reference_steps = true;
    const auto result = accept_seam_feature_targets(chain, targets, options);
    REQUIRE(result.accepted);
    CHECK(seam_feature_continuous(chain, result.positions, result.geometric, options, nullptr, &result.bindings));
    for (size_t i = 0; i < chain.size(); ++i) {
        CHECK(result.locked[i]);
        CHECK_THAT((result.positions[i] - targets[i].target).norm(), Catch::Matchers::WithinAbs(0., 1e-12));
    }
    auto unproved = result.bindings;
    unproved[1].carrier_direction.setZero(); // A geometric flag alone cannot waive the ordinary step bound.
    CHECK_FALSE(seam_feature_continuous(chain, result.positions, result.geometric, options, nullptr, &unproved));
    auto different = result.bindings;
    different.back().edge = 8;
    CHECK_FALSE(seam_feature_continuous(chain, result.positions, result.geometric, options, nullptr, &different));
}

namespace {
std::pair<std::vector<SeamChainPoint>, std::vector<SeamFeatureBinding>> bent_approach()
{
    const Vec2d source[] = {{-7.104284763336182, -38.126434326171875}, {-6.940968036651611, -38.09869384765625},
        {-6.3572869300842285, -37.9473876953125}, {-6.259331226348877, -37.9388542175293},
        {-6.118964195251465, -37.9178581237793}};
    const Vec2d corner[] = {source[0], source[1], {-6.738527774810791, -38.059452056884766},
        {-6.536096096038818, -38.02021026611328}, {-6.3336639404296875, -37.9809684753418}};
    const double weights[] = {1.4831028258, 1.4788977573, .3931916225, .4833962854, .4833977471};
    std::vector<SeamChainPoint> points;
    std::vector<SeamFeatureBinding> targets(5);
    for (size_t i = 0; i < 5; ++i) {
        const Vec2d tangent = i < 2 ? Vec2d::UnitX().eval() : (source[i] - corner[i]).normalized().eval();
        const Vec2d across(-tangent.y(), tangent.x());
        SeamChainPoint p;
        p.width = .42;
        p.anchor_weight = weights[i];
        p.initial = Vec3d(source[i].x(), source[i].y(), .3 + .2 * double(i));
        p.vertices = {{corner[i]}, {corner[i] + 5. * tangent},
            {corner[i] + 5. * tangent + 5. * across}, {corner[i] + 5. * across}};
        p.candidate_arc = p.anchor_arc = (source[i] - corner[i]).norm();
        p.fixed = i < 2; // Existing junction anchors must survive both failed and repaired proposals.
        if (i == 1 || i >= 3) targets[i] = {i == 1 ? SeamFeatureState::Locked : SeamFeatureState::Proposed,
            9, 0, Vec3d(corner[i].x(), corner[i].y(), p.initial.z()), 0.};
        points.push_back(p);
    }
    return {points, targets};
}
} // namespace

TEST_CASE("A turn repair retains an edge approach without moving its existing anchors", "[SeamPlacer][SeamFeature][Regression]")
{
    const auto fixture = bent_approach();
    const auto &points = fixture.first;
    const auto &targets = fixture.second;
    SeamFeatureOptions options;
    options.bridge_repair_turns = false;
    CHECK_FALSE(accept_seam_feature_targets(points, targets, options).accepted);
    options.bridge_repair_turns = true;
    const auto result = accept_seam_feature_targets(points, targets, options);
    REQUIRE(result.accepted);
    CHECK(seam_feature_continuous(points, result.positions, result.geometric, options));
    for (size_t i : {size_t(0), size_t(1)})
        CHECK_THAT((result.positions[i] - points[i].initial).norm(), Catch::Matchers::WithinAbs(0., 1e-12));
    for (size_t i : {size_t(3), size_t(4)}) {
        CHECK(result.locked[i]);
        CHECK_THAT((result.positions[i] - targets[i].target).norm(), Catch::Matchers::WithinAbs(0., 1e-12));
    }
    CHECK_FALSE(result.locked[2]);
    CHECK((result.positions[2] - points[2].initial).norm() <= options.bridge_motion_widths * points[2].width + 1e-5);
}

TEST_CASE("A turn repair cannot cross a painted boundary to keep a proposed approach", "[SeamPlacer][SeamFeature][Regression]")
{
    auto fixture = bent_approach();
    auto &points = fixture.first;
    points[2].lower_arc = points[2].candidate_arc - .1; // The useful correction is outside this enforced interval.
    points[2].upper_arc = points[2].candidate_arc;
    const auto result = accept_seam_feature_targets(points, fixture.second);
    CHECK_FALSE(result.accepted);
    for (size_t i = 0; i < points.size(); ++i)
        CHECK_THAT((result.positions[i] - points[i].initial).norm(), Catch::Matchers::WithinAbs(0., 1e-12));
}
