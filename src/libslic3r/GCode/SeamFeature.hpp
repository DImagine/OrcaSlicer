#pragma once

#include "SeamChain.hpp"
#include "libslic3r/Point.hpp"
#include <admesh/stl.h>
#include <cstddef>
#include <vector>

namespace Slic3r::SeamPlacerImpl {

// Experimental controls are independent of original candidate ranking and the 55-degree snap.
struct SeamFeatureOptions {
    double corridor_widths = 4.;
    double minimum_angle = 15.;
    double minimum_edge_widths = 1.;
    double trace_angle = 3.;
    double offset_min_widths = .25;
    double offset_max_widths = 1.5;
    double offset_difference_widths = .12;
    double geometry_tolerance = .012;
    double node_widths = 2.;
    double horizontal_normal_min = .15;
    double capture_widths = .75;
    size_t evidence_layers = 3;
    size_t evidence_window = 5;
    double evidence_span_widths = .5;
    double route_margin = .05;
    double unbound_cost = 1.;
    double transition_span = 4.;
    double transition_bending = .05;
    double bridge_motion_widths = .75;
    double arc_chart_fraction = .49;
    double step_growth = .25;
    double step_allowance_widths = .15;
    double turn_allowance = 20.;
    double support_allowance_widths = .1;
    bool carrier_reference_steps = true;
    bool bridge_repair_turns = true;
    size_t bridge_repair_samples = 32;
    size_t max_iterations = 128;
    size_t bridge_repair_sweeps = 8;
    size_t max_chain_points = 1500;
    size_t max_edge_queries = 1000000;
    size_t max_candidates = 32;
    double solve_tolerance = 1e-7;
};

enum class SeamFeatureState { None, Ambiguous, InsufficientEvidence, Proposed, Locked, TransitionRejected, Budget };
struct SeamFeatureBinding {
    SeamFeatureState state = SeamFeatureState::None;
    size_t edge = size_t(-1);
    size_t vertex = size_t(-1);
    Vec3d target = Vec3d::Zero();
    double arc = 0.;
    Vec3d carrier_direction = Vec3d::Zero();
};
struct SeamFeatureResult {
    bool accepted = false;
    std::vector<Vec3d> positions;
    std::vector<bool> locked;
    std::vector<bool> geometric;
    std::vector<SeamFeatureBinding> bindings;
};

// Collect geometry for every chain before proposing any immutable binding.
std::vector<SeamFeatureResult> bind_seam_features(const indexed_triangle_set &mesh,
    const std::vector<std::vector<SeamChainPoint>> &chains, const SeamFeatureOptions &options = {});
// Transition construction is transactional: infeasible proposed runs remain at their original positions.
SeamFeatureResult accept_seam_feature_targets(const std::vector<SeamChainPoint> &points,
    std::vector<SeamFeatureBinding> bindings, const SeamFeatureOptions &options = {});
// Recheck complete later optimizer results against the same original line and frozen binding mask.
bool seam_feature_continuous(const std::vector<SeamChainPoint> &original, const std::vector<Vec3d> &positions,
    const std::vector<bool> &geometric, const SeamFeatureOptions &options = {}, size_t *failed = nullptr,
    const std::vector<SeamFeatureBinding> *bindings = nullptr);

} // namespace Slic3r::SeamPlacerImpl
