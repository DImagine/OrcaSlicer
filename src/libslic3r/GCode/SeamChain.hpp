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
    // An excluded inclination freezes motion without becoming a geometric target.
    bool projection_hold = false;
    double projection_angle = std::numeric_limits<double>::infinity();
    // Frozen model-surface patch from the preceding length-only pass; zero confidence disables its preference.
    Vec3d surface_normal = Vec3d::Zero();
    double surface_confidence = 0.;
    // A sharp boundary may shrink the chart below the nominal projection radius.
    double surface_radius = std::numeric_limits<double>::infinity();
};

enum class SeamChainStatus { Disabled, Invalid, NoFreedom, Unchanged, Improved, Budget, Failed };

// This response belongs only to post-fit optimization, never to original candidate selection or fitting.
enum class SeamChainAngleCurve { Legacy, Smoothstep, Linear, Cubic };
struct SeamChainAngleOptions {
    SeamChainAngleCurve curve = SeamChainAngleCurve::Legacy;
    double minimum_degrees = 15.;
    double full_degrees = 30.;
    double painted_center_weight = 3.;
};

// Preserve signed angle preference in base_weight; apply one magnitude response to both paint types.
double seam_chain_anchor_weight(double base_weight, double legacy_factor, double angle_degrees,
                               bool painted, const SeamChainAngleOptions &options = {});

struct SeamChainOptions {
    double length_weight = 0.15;
    double vertical_weight = 0.;
    double projection_span = 4.;
    // Zero uses an arc-length window; three uses the immediate neighboring layers.
    size_t direction_layers = 0;
    bool geometric_probe = false;
    // Projection settings are centralized here; the ordinary length pass remains unchanged.
    double relative_length_growth = 0.015;
    double max_projection_angle = 10.;
    double max_face_tilt = 45.;
    double sharp_surface_angle = 55.;
    double max_motion_widths = 4.;
    double boundary_taper_fraction = 0.25;
    bool projection_gate = false;
    bool hold_free_center = false;
    size_t segments = 1;
    size_t max_iterations = 100;
    size_t max_evaluations = 5000;
    static SeamChainOptions projected_defaults() {
        SeamChainOptions options;
        options.vertical_weight = 1.; // Keep the accepted vertical multiplier independent of the length weight.
        options.direction_layers = 3;
        options.geometric_probe = options.projection_gate = options.hold_free_center = true;
        return options;
    }
};

struct SeamChainResult {
    SeamChainStatus status = SeamChainStatus::Invalid;
    std::vector<Vec3d> positions;
    std::vector<double> arcs;
    double initial_energy = 0.;
    double energy = 0.;
    double initial_length = 0.;
    double length = 0.;
    // Keep the objective components separate so parameter trials can explain accepted tradeoffs.
    double preference_cost = 0.;
    double length_cost = 0.;
    double projection_cost = 0.;
    size_t iterations = 0;
    // Diagnostics never replace the committed energy/positions when the chain falls back.
    double trial_energy = 0.;
    double trial_length = 0.;
    size_t length_rejections = 0;
    double stationarity = std::numeric_limits<double>::infinity();
    bool resolution_converged = false;
    int qp_failure = 0;
};

// Only Improved commits positions. Every other status returns the complete original chain.
SeamChainResult optimize_seam_chain(const std::vector<SeamChainPoint> &points, const SeamChainOptions &options = {});
// Freeze eligibility and the optional free-chain midpoint once, after the ordinary length pass.
void prepare_seam_projection(std::vector<SeamChainPoint> &points, const SeamChainOptions &options);

} // namespace Slic3r::SeamPlacerImpl
