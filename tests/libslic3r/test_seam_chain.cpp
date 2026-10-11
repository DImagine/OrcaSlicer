#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/GCode/SeamChain.hpp"
#include "libslic3r/libslic3r.h"
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

using namespace Slic3r;

TEST_CASE("Length optimization shortens a smooth seam while keeping its fixed ends", "[SeamPlacer][SeamChain][Regression]")
{
    using namespace SeamPlacerImpl;
    std::vector<SeamChainPoint> chain;
    for (size_t i = 0; i <= 40; ++i) {
        SeamChainPoint point;
        const double z = double(i) * 0.5;
        const double x = 0.6 * std::sin(2. * PI * z / 20.);
        point.vertices = {{Vec2d(-10, 0)}, {Vec2d(10, 0)}, {Vec2d(10, 10)}, {Vec2d(-10, 10)}};
        point.initial = Vec3d(x, 0., z);
        point.candidate_arc = point.anchor_arc = x + 10.;
        point.fixed = i == 0 || i == 20 || i == 40;
        chain.push_back(std::move(point));
    }
    SeamChainOptions options;
    options.segments = 5;
    const auto result = optimize_seam_chain(chain, options);
    REQUIRE(result.status == SeamChainStatus::Improved);
    CHECK(result.length < result.initial_length - 0.05);
    CHECK(result.energy < result.initial_energy);
    CHECK((result.positions.front() - chain.front().initial).norm() < 1e-12);
    CHECK((result.positions.back() - chain.back().initial).norm() < 1e-12);
    CHECK((result.positions[20] - chain[20].initial).norm() < 1e-12);
    for (const auto &position : result.positions) CHECK(std::abs(position.y()) < 1e-10);
}

TEST_CASE("Length optimization respects narrow painted ranges on every layer", "[SeamPlacer][SeamChain][Regression]")
{
    using namespace SeamPlacerImpl;
    std::vector<SeamChainPoint> chain;
    for (size_t i = 0; i <= 30; ++i) {
        SeamChainPoint point;
        const double z = double(i) * 0.5;
        const double x = 0.15 * std::sin(2. * PI * z / 15.);
        point.vertices = {{Vec2d(-10, 0)}, {Vec2d(10, 0)}, {Vec2d(10, 10)}, {Vec2d(-10, 10)}};
        point.initial = Vec3d(x, 0., z);
        point.candidate_arc = point.anchor_arc = x + 10.;
        // The feasible tube changes along the height; some layers force the path to remain off center.
        point.lower_arc = i >= 12 && i <= 18 ? 10.03 : 9.8;
        point.upper_arc = 10.2;
        if (point.initial.x() < point.lower_arc - 10.) point.initial.x() = point.lower_arc - 10.;
        chain.push_back(std::move(point));
    }
    SeamChainOptions options;
    options.segments = 5;
    const auto result = optimize_seam_chain(chain, options);
    REQUIRE(result.status == SeamChainStatus::Improved);
    for (size_t i = 0; i < chain.size(); ++i) {
        CHECK(result.arcs[i] >= chain[i].lower_arc - 1e-8);
        CHECK(result.arcs[i] <= chain[i].upper_arc + 1e-8);
        CHECK(std::abs(result.positions[i].y()) < 1e-10);
    }
}

TEST_CASE("An exhausted seam optimization budget returns every original point", "[SeamPlacer][SeamChain][Regression]")
{
    using namespace SeamPlacerImpl;
    const size_t budget = GENERATE(size_t(0), size_t(1));
    std::vector<SeamChainPoint> chain;
    for (size_t i = 0; i <= 20; ++i) {
        SeamChainPoint point;
        const double x = 0.6 * std::sin(double(i) * 0.3);
        point.vertices = {{Vec2d(-10, 0)}, {Vec2d(10, 0)}, {Vec2d(10, 10)}, {Vec2d(-10, 10)}};
        point.initial = Vec3d(x, 0., double(i) * 0.5);
        point.candidate_arc = point.anchor_arc = x + 10.;
        chain.push_back(std::move(point));
    }
    SeamChainOptions options;
    options.segments = 4;
    options.max_iterations = budget;
    const auto result = optimize_seam_chain(chain, options);
    CHECK(result.status == SeamChainStatus::Budget);
    for (size_t i = 0; i < chain.size(); ++i) CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
}

TEST_CASE("Disabled and fully anchored seam chains keep the original positions", "[SeamPlacer][SeamChain]")
{
    using namespace SeamPlacerImpl;
    const bool disabled = GENERATE(false, true);
    std::vector<SeamChainPoint> chain;
    for (size_t i = 0; i < 12; ++i) {
        SeamChainPoint point;
        point.vertices = {{Vec2d(-10, 0)}, {Vec2d(10, 0)}, {Vec2d(10, 10)}, {Vec2d(-10, 10)}};
        point.initial = Vec3d(double(i) * 0.02, 0., double(i) * 0.2);
        point.fixed = true;
        chain.push_back(std::move(point));
    }
    SeamChainOptions options;
    options.length_weight = disabled ? 0. : 0.1;
    const auto result = optimize_seam_chain(chain, options);
    CHECK(result.status == (disabled ? SeamChainStatus::Disabled : SeamChainStatus::NoFreedom));
    for (size_t i = 0; i < chain.size(); ++i) CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
}

TEST_CASE("Length optimization crosses polygon vertices with a painted end and a free end", "[SeamPlacer][SeamChain][Regression]")
{
    using namespace SeamPlacerImpl;
    std::vector<SeamChainPoint> chain;
    constexpr size_t sides = 48;
    constexpr double radius = 10.;
    const double edge_length = 2. * radius * std::sin(PI / double(sides));
    for (size_t i = 0; i <= 60; ++i) {
        SeamChainPoint point;
        for (size_t j = 0; j < sides; ++j) {
            const double angle = 2. * PI * double(j) / double(sides);
            point.vertices.push_back({Vec2d(radius * std::cos(angle), radius * std::sin(angle))});
        }
        const double displacement = 0.2 * std::sin(2. * PI * double(i) / 60.);
        // The original path crosses the contour's coordinate cut and actual polygon vertices.
        const size_t edge = displacement < 0. ? sides - 1 : 0;
        const double fraction = displacement < 0. ? 1. + displacement / edge_length : displacement / edge_length;
        const Vec2d xy = point.vertices[edge].position + fraction *
            (point.vertices[(edge + 1) % sides].position - point.vertices[edge].position);
        point.initial = Vec3d(xy.x(), xy.y(), 0.5 * double(i));
        point.candidate_arc = displacement;
        point.anchor_arc = 0.;
        if (i < 20) {
            point.anchor_weight = 3.;
            point.lower_arc = -0.25;
            point.upper_arc = 0.25;
        }
        chain.push_back(std::move(point));
    }
    SeamChainOptions options;
    options.segments = 8;
    const auto result = optimize_seam_chain(chain, options);
    INFO("iterations=" << result.iterations << " residual=" << result.stationarity << " qp_failure=" << result.qp_failure
         << " initial=" << result.initial_energy << " trial=" << result.trial_energy);
    REQUIRE(result.status == SeamChainStatus::Improved);
    CHECK(result.length < result.initial_length - 0.002);
    for (size_t i = 0; i < chain.size(); ++i) {
        CHECK(result.arcs[i] >= chain[i].lower_arc - 1e-8);
        CHECK(result.arcs[i] <= chain[i].upper_arc + 1e-8);
    }
}

TEST_CASE("Evaluation exhaustion and invalid geometry preserve the entire original seam", "[SeamPlacer][SeamChain][Regression]")
{
    using namespace SeamPlacerImpl;
    const bool invalid = GENERATE(false, true);
    std::vector<SeamChainPoint> chain;
    for (size_t i = 0; i < 15; ++i) {
        SeamChainPoint point;
        // A duplicate vertex verifies that zero-length contour edges do not break parameterization.
        point.vertices = {{Vec2d(-10, 0)}, {Vec2d(-10, 0)}, {Vec2d(10, 0)}, {Vec2d(10, 10)}, {Vec2d(-10, 10)}};
        const double x = 0.4 * std::sin(double(i) * 0.3);
        point.initial = Vec3d(x, 0., double(i) * 0.25);
        point.candidate_arc = point.anchor_arc = x + 10.;
        chain.push_back(std::move(point));
    }
    SeamChainOptions options;
    options.segments = 3;
    if (invalid) chain[7].initial.y() = -0.01; // An initial point outside its contour must never be silently projected.
    else options.max_evaluations = GENERATE(size_t(0), size_t(1));
    const auto result = optimize_seam_chain(chain, options);
    CHECK(result.status == (invalid ? SeamChainStatus::Invalid : SeamChainStatus::Budget));
    for (size_t i = 0; i < chain.size(); ++i) CHECK((result.positions[i] - chain[i].initial).norm() < 1e-12);
}

