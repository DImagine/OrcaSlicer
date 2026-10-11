#pragma once

#include "libslic3r/Point.hpp"
#include <cstddef>
#include <limits>
#include <vector>

namespace Slic3r::SeamPlacerImpl {

// Samples keep the old candidate scores; geometry is interpolated on the closed layer contour.
struct SeamChainVertex {
    Vec2d position;
    double quality = 0.;
    double overhang = 0.;
};

struct SeamChainPoint {
    std::vector<SeamChainVertex> vertices;
    Vec3d initial = Vec3d::Zero();
    double candidate_arc = 0.;
    double anchor_arc = 0.;
    double anchor_weight = 0.;
    double width = 0.4;
    double lower_arc = -std::numeric_limits<double>::infinity();
    double upper_arc = std::numeric_limits<double>::infinity();
    bool fixed = false;
};

enum class SeamChainStatus { Disabled, Invalid, NoFreedom, Unchanged, Improved, Budget, Failed };

struct SeamChainOptions {
    double length_weight = 0.3;
    size_t segments = 1;
    size_t max_iterations = 100;
    size_t max_evaluations = 5000;
};

struct SeamChainResult {
    SeamChainStatus status = SeamChainStatus::Invalid;
    std::vector<Vec3d> positions;
    std::vector<double> arcs;
    double initial_energy = 0.;
    double energy = 0.;
    double initial_length = 0.;
    double length = 0.;
    size_t iterations = 0;
    // Diagnostics never replace the committed energy/positions when the chain falls back.
    double trial_energy = 0.;
    double stationarity = std::numeric_limits<double>::infinity();
    bool resolution_converged = false;
    int qp_failure = 0;
};

// Only Improved commits positions. Every other status returns the complete original chain.
SeamChainResult optimize_seam_chain(const std::vector<SeamChainPoint> &points, const SeamChainOptions &options = {});

} // namespace Slic3r::SeamPlacerImpl
