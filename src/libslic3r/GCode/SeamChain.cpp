#include "SeamChain.hpp"
#include "libslic3r/Geometry/Bicubic.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <initializer_list>
#include <numeric>
#include <utility>
#include <vector>

namespace Slic3r::SeamPlacerImpl {
namespace {
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
constexpr double geometry_tolerance = 1e-5;

struct Contour {
    const SeamChainPoint *point = nullptr;
    std::vector<double> arc;
    double total = 0.;
    double initial_arc = 0.;
    double anchor_quality = 0.;
    double anchor_overhang = 0.;
};

struct Sample {
    Vec3d position;
    Vec3d tangent;
    double quality;
    double quality_slope;
    double overhang;
    double overhang_slope;
};

Sample sample(const Contour &contour, double arc)
{
    const auto &vertices = contour.point->vertices;
    arc = std::fmod(arc, contour.total);
    if (arc < 0.) arc += contour.total;
    // upper_bound skips repeated vertices and gives the outgoing nonzero edge at a vertex.
    size_t edge = size_t(std::upper_bound(contour.arc.begin(), contour.arc.end(), arc) - contour.arc.begin() - 1);
    edge = std::min(edge, vertices.size() - 1);
    const auto &a = vertices[edge];
    const auto &b = vertices[(edge + 1) % vertices.size()];
    const double length = contour.arc[edge + 1] - contour.arc[edge];
    const double t = (arc - contour.arc[edge]) / length;
    const Vec2d xy = a.position + t * (b.position - a.position);
    Vec2d tangent = (b.position - a.position) / length;
    double q_slope = (b.quality - a.quality) / length;
    double o_slope = (b.overhang - a.overhang) / length;
    if (arc - contour.arc[edge] <= 1e-7) {
        // At a mesh vertex use both one-sided models; the objective still samples the exact polyline.
        size_t previous = (edge + vertices.size() - 1) % vertices.size();
        for (size_t visited = 0; visited < vertices.size(); ++visited) {
            const double previous_length = contour.arc[previous + 1] - contour.arc[previous];
            if (previous_length > 0.) {
                tangent = 0.5 * (tangent + (a.position - vertices[previous].position) / previous_length);
                q_slope = 0.5 * (q_slope + (a.quality - vertices[previous].quality) / previous_length);
                o_slope = 0.5 * (o_slope + (a.overhang - vertices[previous].overhang) / previous_length);
                break;
            }
            previous = (previous + vertices.size() - 1) % vertices.size();
        }
    }
    return {Vec3d(xy.x(), xy.y(), contour.point->initial.z()), Vec3d(tangent.x(), tangent.y(), 0.),
            a.quality + t * (b.quality - a.quality), q_slope,
            a.overhang + t * (b.overhang - a.overhang), o_slope};
}

bool prepare(const SeamChainPoint &point, Contour &contour)
{
    if (point.vertices.size() < 3 || !point.initial.allFinite() || !(point.width > 0.) ||
        !std::isfinite(point.width) || !std::isfinite(point.anchor_weight) || point.anchor_weight < 0. ||
        !std::isfinite(point.anchor_arc) || !std::isfinite(point.candidate_arc) ||
        std::isnan(point.lower_arc) || std::isnan(point.upper_arc) || point.lower_arc > point.upper_arc)
        return false;
    contour.point = &point;
    contour.arc.push_back(0.);
    double nearest = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < point.vertices.size(); ++i) {
        const auto &a = point.vertices[i];
        const auto &b = point.vertices[(i + 1) % point.vertices.size()];
        if (!a.position.allFinite() || !std::isfinite(a.quality) || !std::isfinite(a.overhang)) return false;
        const Vec2d edge = b.position - a.position;
        const double length = edge.norm();
        if (length > 0.) {
            const double t = std::clamp((point.initial.head<2>() - a.position).dot(edge) / edge.squaredNorm(), 0., 1.);
            const double distance = (point.initial.head<2>() - a.position - t * edge).norm();
            if (distance < nearest) {
                nearest = distance;
                contour.initial_arc = contour.arc.back() + t * length;
            }
        }
        contour.arc.push_back(contour.arc.back() + length);
    }
    contour.total = contour.arc.back();
    if (!(contour.total > 0.) || nearest > geometry_tolerance) return false;
    // Choose the unwrapped copy of the initial point belonging to its permitted painted/blocked arc.
    contour.initial_arc += std::round((point.candidate_arc - contour.initial_arc) / contour.total) * contour.total;
    if (contour.initial_arc < point.lower_arc - geometry_tolerance ||
        contour.initial_arc > point.upper_arc + geometry_tolerance) return false;
    const Sample anchor = sample(contour, point.anchor_arc);
    contour.anchor_quality = anchor.quality;
    contour.anchor_overhang = anchor.overhang;
    return true;
}

struct Evaluation {
    double energy = 0.;
    double length = 0.;
    Vector gradient;
    Matrix hessian;
    std::vector<Vec3d> positions;
    Vector arcs;
};

Evaluation evaluate(const std::vector<Contour> &contours, const Matrix &basis, const Vector &variables,
                    double width, double lambda, bool derivatives)
{
    const size_t count = contours.size();
    const Eigen::Index dimensions = variables.size();
    Evaluation value;
    value.gradient = Vector::Zero(dimensions);
    value.hessian = Matrix::Zero(dimensions, dimensions);
    value.arcs = width * basis * variables;
    value.positions.reserve(count);
    Matrix last_jacobian;
    for (size_t i = 0; i < count; ++i) {
        const Contour &contour = contours[i];
        const SeamChainPoint &point = *contour.point;
        value.arcs[Eigen::Index(i)] += contour.initial_arc;
        const Sample current = sample(contour, value.arcs[Eigen::Index(i)]);
        // Fixed points retain the exact original coordinates rather than their numerically reconstructed foot.
        value.positions.push_back(point.fixed ? point.initial : current.position);
        double height_weight = 0.;
        if (i > 0) height_weight += 0.5 * (point.initial.z() - contours[i - 1].point->initial.z());
        if (i + 1 < count) height_weight += 0.5 * (contours[i + 1].point->initial.z() - point.initial.z());
        height_weight /= width;
        // Constant costs at frozen rows must not inflate relative convergence tolerances.
        if (point.fixed || basis.row(Eigen::Index(i)).norm() < 1e-14) height_weight = 0.;
        double distance = value.arcs[Eigen::Index(i)] - point.anchor_arc;
        // Periodic free contours have no artificial movement radius or privileged coordinate cut.
        if (!std::isfinite(point.lower_arc) && !std::isfinite(point.upper_arc))
            distance = std::remainder(distance, contour.total);
        const double q_residual = std::max(0., (current.quality - contour.anchor_quality - 0.3) / 0.3);
        const double o_tolerance = 0.1 * point.width;
        const double o_residual = std::max(0., (current.overhang - contour.anchor_overhang - o_tolerance) / o_tolerance);
        const double anchor_factor = point.anchor_weight / (point.width * point.width);
        value.energy += height_weight * (anchor_factor * distance * distance + q_residual * q_residual + o_residual * o_residual);
        const Vector row = width * basis.row(Eigen::Index(i)).transpose();
        const Matrix jacobian = current.tangent * row.transpose();
        if (derivatives) {
            const double q_slope = q_residual > 0. ? current.quality_slope / 0.3 : 0.;
            const double o_slope = o_residual > 0. ? current.overhang_slope / o_tolerance : 0.;
            value.gradient += 2. * height_weight * (anchor_factor * distance + q_residual * q_slope + o_residual * o_slope) * row;
            value.hessian += 2. * height_weight * (anchor_factor + q_slope * q_slope + o_slope * o_slope) * row * row.transpose();
        }
        if (i > 0) {
            const Vec3d delta = value.positions[i] - value.positions[i - 1];
            const double length = delta.norm();
            value.length += length;
            value.energy += lambda / width * delta.head<2>().squaredNorm() / (length + delta.z());
            if (derivatives) {
                const Matrix difference = jacobian - last_jacobian;
                value.gradient += lambda / (width * length) * difference.transpose() * delta;
                // A positive quadratic majorizer of length avoids an indefinite Newton system.
                value.hessian += lambda / (width * length) * difference.transpose() * difference;
            }
        }
        last_jacobian = jacobian;
    }
    return value;
}

// Convex QP with a feasible zero start. Every working-set loop has a hard iteration budget.
bool quadratic_step(const Matrix &hessian, const Vector &gradient, const Matrix &constraints,
                    const Vector &limits, Vector &step, int *failure = nullptr)
{
    auto fail = [&](int reason) { if (failure) *failure = reason; return false; };
    const Eigen::Index dimensions = gradient.size();
    step = Vector::Zero(dimensions);
    std::vector<Eigen::Index> active;
    for (size_t iteration = 0; iteration < 200; ++iteration) {
        Matrix active_rows(Eigen::Index(active.size()), dimensions);
        for (size_t j = 0; j < active.size(); ++j) active_rows.row(Eigen::Index(j)) = constraints.row(active[j]);
        Matrix tangent_space = Matrix::Identity(dimensions, dimensions);
        if (!active.empty()) {
            // Solve in the constraint null space: a KKT factorization squares the conditioning of nearby rows.
            const Eigen::FullPivLU<Matrix> factor(active_rows);
            if (factor.dimensionOfKernel() == 0) tangent_space.resize(dimensions, 0);
            else tangent_space = factor.kernel();
        }
        const Vector force = gradient + hessian * step;
        Vector direction = Vector::Zero(dimensions);
        if (tangent_space.cols() > 0) {
            const Matrix reduced = tangent_space.transpose() * hessian * tangent_space;
            const Eigen::LDLT<Matrix> factor(reduced);
            if (factor.info() != Eigen::Success || !factor.isPositive()) return fail(1);
            direction = tangent_space * factor.solve(-tangent_space.transpose() * force);
        }
        if (!direction.allFinite()) return fail(2);
        if (direction.lpNorm<Eigen::Infinity>() < 1e-10) {
            Vector multipliers = Vector::Zero(Eigen::Index(active.size()));
            if (!active.empty()) multipliers = active_rows.transpose().colPivHouseholderQr().solve(-force);
            if (!multipliers.allFinite() || (active_rows.transpose() * multipliers + force).norm() > 1e-7 * std::max(1., force.norm()))
                return fail(2);
            // Nonnegative multipliers certify optimality for the current convex subproblem.
            Eigen::Index remove = -1;
            double minimum = -1e-8;
            for (size_t j = 0; j < active.size(); ++j)
                if (multipliers[Eigen::Index(j)] < minimum) {
                    minimum = multipliers[Eigen::Index(j)];
                    remove = Eigen::Index(j);
                }
            if (remove < 0) return (constraints * step - limits).maxCoeff() <= 1e-7 ? true : fail(3);
            active.erase(active.begin() + remove);
            continue;
        }
        double fraction = 1.;
        Eigen::Index blocking = -1;
        for (Eigen::Index row = 0; row < constraints.rows(); ++row) {
            if (std::find(active.begin(), active.end(), row) != active.end()) continue;
            const double derivative = constraints.row(row).dot(direction);
            if (derivative <= 1e-10 * std::max(1., constraints.row(row).norm() * direction.norm())) continue;
            const double allowed = std::max(0., (limits[row] - constraints.row(row).dot(step)) / derivative);
            if (allowed < fraction) { fraction = allowed; blocking = row; }
        }
        step += fraction * direction;
        if (blocking >= 0) active.push_back(blocking);
    }
    return fail(4);
}
} // namespace

SeamChainResult optimize_seam_chain(const std::vector<SeamChainPoint> &points, const SeamChainOptions &options)
{
    SeamChainResult result;
    for (size_t i = 0; i < points.size(); ++i) {
        result.positions.push_back(points[i].initial);
        if (i > 0) result.initial_length += (points[i].initial - points[i - 1].initial).norm();
    }
    result.length = result.initial_length;
    if (options.length_weight == 0.) { result.status = SeamChainStatus::Disabled; return result; }
    if (points.size() < 2 || options.segments == 0 || !std::isfinite(options.length_weight) || options.length_weight < 0.) return result;
    std::vector<Contour> contours(points.size());
    std::vector<double> widths;
    std::vector<size_t> fixed;
    for (size_t i = 0; i < points.size(); ++i) {
        if (!prepare(points[i], contours[i]) || (i > 0 && points[i].initial.z() <= points[i - 1].initial.z())) return result;
        widths.push_back(points[i].width);
        result.arcs.push_back(contours[i].initial_arc);
        if (points[i].fixed) fixed.push_back(i);
    }
    std::sort(widths.begin(), widths.end());
    const double width = widths[widths.size() / 2];
    const size_t segments = std::min(options.segments, points.size() - 1);
    const double segment_height = (points.back().initial.z() - points.front().initial.z()) / double(segments);
    Matrix basis = Matrix::Zero(Eigen::Index(points.size()), Eigen::Index(segments + 1));
    for (size_t i = 0; i < points.size(); ++i) {
        const double z = (points[i].initial.z() - points.front().initial.z()) / segment_height;
        const int first = int(std::floor(z)) - 1;
        for (int j = first; j < first + 4; ++j)
            basis(Eigen::Index(i), std::clamp(j, 0, int(segments))) += Geometry::CubicBSplineKernel<double>::kernel(double(j) - z);
    }
    if (!fixed.empty()) {
        Matrix fixed_basis(Eigen::Index(fixed.size()), basis.cols());
        for (size_t j = 0; j < fixed.size(); ++j) fixed_basis.row(Eigen::Index(j)) = basis.row(Eigen::Index(fixed[j]));
        Eigen::FullPivLU<Matrix> factor(fixed_basis);
        if (factor.dimensionOfKernel() == 0) { result.status = SeamChainStatus::NoFreedom; return result; }
        const Matrix null_space = factor.kernel();
        basis = (basis * null_space).eval();
        // Eliminate roundoff in fixed rows so no objective evaluation can move an anchored seam.
        for (size_t i : fixed) basis.row(Eigen::Index(i)).setZero();
    }
    if (basis.norm() < 1e-12) { result.status = SeamChainStatus::NoFreedom; return result; }
    const Eigen::Index dimensions = basis.cols();
    std::vector<Vector> poll_targets;
    for (Eigen::Index column = 0; column < dimensions; ++column) {
        Vector target = Vector::Zero(dimensions);
        target[column] = 1.;
        poll_targets.push_back(std::move(target));
    }
    // Coherent modes can cross small facet barriers that trap a search moving one local coefficient.
    const Eigen::LDLT<Matrix> basis_fit(basis.transpose() * basis);
    auto add_mode = [&](const Vector &profile) {
        Vector target = basis_fit.solve(basis.transpose() * profile);
        if (!target.allFinite()) return;
        const double amplitude = (basis * target).lpNorm<Eigen::Infinity>();
        if (amplitude > 1e-12) poll_targets.push_back(target / amplitude);
    };
    add_mode(Vector::Ones(basis.rows()));
    for (int mode = 1; mode <= 3; ++mode) {
        Vector profile(basis.rows());
        for (size_t i = 0; i < points.size(); ++i)
            profile[Eigen::Index(i)] = std::cos(double(mode) * PI *
                (points[i].initial.z() - points.front().initial.z()) /
                (points.back().initial.z() - points.front().initial.z()));
        add_mode(profile);
    }
    for (size_t first = 0; first < points.size();) {
        const auto free_contour = [&](size_t i) { return !points[i].fixed &&
            !std::isfinite(points[i].lower_arc) && !std::isfinite(points[i].upper_arc); };
        if (!free_contour(first)) { ++first; continue; }
        size_t end = first + 1;
        while (end < points.size() && free_contour(end)) ++end;
        if (end - first >= 6) {
            Vector profile = Vector::Zero(basis.rows());
            profile.segment(Eigen::Index(first), Eigen::Index(end - first)).setOnes();
            add_mode(profile);
        }
        first = end;
    }
    std::vector<Vector> rows;
    std::vector<double> upper;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].fixed) continue;
        if (std::isfinite(points[i].upper_arc)) {
            rows.push_back(basis.row(Eigen::Index(i)).transpose());
            upper.push_back(std::max(0., (points[i].upper_arc - contours[i].initial_arc) / width));
        }
        if (std::isfinite(points[i].lower_arc)) {
            rows.push_back(-basis.row(Eigen::Index(i)).transpose());
            upper.push_back(std::max(0., (contours[i].initial_arc - points[i].lower_arc) / width));
        }
    }
    Matrix constraints(Eigen::Index(rows.size()), dimensions);
    Vector limits(Eigen::Index(rows.size()));
    for (size_t i = 0; i < rows.size(); ++i) { constraints.row(Eigen::Index(i)) = rows[i].transpose(); limits[Eigen::Index(i)] = upper[i]; }
    Vector variables = Vector::Zero(dimensions);
    if (options.max_evaluations == 0) { result.status = SeamChainStatus::Budget; return result; }
    Evaluation initial = evaluate(contours, basis, variables, width, options.length_weight, true);
    Evaluation current = initial;
    result.initial_energy = result.energy = initial.energy;
    if (!std::isfinite(initial.energy)) return result;
    size_t evaluations = 1;
    double damping = 1e-3;
    bool converged = false;
    for (size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        result.iterations = iteration + 1;
        if (evaluations >= options.max_evaluations) break;
        bool accepted = false;
        bool stationary_model = false;
        bool small_progress = false;
        if (!current.gradient.allFinite() || !current.hessian.allFinite()) { result.status = SeamChainStatus::Failed; return result; }
        for (size_t attempt = 0; attempt < 8; ++attempt) {
            Matrix hessian = current.hessian;
            for (Eigen::Index j = 0; j < dimensions; ++j) hessian(j, j) += damping * std::max(1., hessian(j, j));
            // Bounds on B*step limit one deformation, not the total movement along the layer contour.
            Matrix step_constraints(constraints.rows() + 2 * basis.rows(), dimensions);
            Vector step_limits(step_constraints.rows());
            step_constraints.topRows(constraints.rows()) = constraints;
            step_limits.head(constraints.rows()) = limits - constraints * variables;
            step_constraints.middleRows(constraints.rows(), basis.rows()) = basis;
            step_constraints.bottomRows(basis.rows()) = -basis;
            step_limits.tail(2 * basis.rows()).setConstant(0.5);
            Vector step;
            if (!quadratic_step(hessian, current.gradient, step_constraints, step_limits, step)) {
                damping *= 10.;
                continue;
            }
            const double movement = width * (basis * step).lpNorm<Eigen::Infinity>();
            // Inspect an undamped projected gradient separately: damping alone cannot certify convergence.
            Vector projected;
            const bool projected_ok = quadratic_step(Matrix::Identity(dimensions, dimensions), current.gradient,
                                                      step_constraints, step_limits, projected);
            const double stationarity = projected_ok ? projected.norm() / std::max(1., initial.energy) : std::numeric_limits<double>::infinity();
            result.stationarity = stationarity;
            if (stationarity <= 1e-6 && movement <= 0.001) { stationary_model = true; break; }
            for (size_t backtrack = 0; backtrack < 12; ++backtrack) {
                if (++evaluations > options.max_evaluations) break;
                const Vector trial_variables = variables + std::ldexp(1., -int(backtrack)) * step;
                if (constraints.rows() > 0 && (constraints * trial_variables - limits).maxCoeff() > 1e-8) continue;
                Evaluation trial = evaluate(contours, basis, trial_variables, width, options.length_weight, false);
                if (!std::isfinite(trial.energy) || trial.energy >= current.energy - 1e-14 * std::max(1., current.energy)) continue;
                small_progress = current.energy - trial.energy <= 1e-8 * std::max(1., initial.energy) &&
                                 std::ldexp(movement, -int(backtrack)) <= 0.001;
                variables = trial_variables;
                if (evaluations >= options.max_evaluations) { result.status = SeamChainStatus::Budget; return result; }
                current = evaluate(contours, basis, variables, width, options.length_weight, true);
                result.trial_energy = current.energy;
                ++evaluations;
                damping = std::max(1e-8, damping * 0.5);
                accepted = true;
                break;
            }
            if (accepted || evaluations >= options.max_evaluations) break;
            damping *= 10.;
        }
        if ((!accepted || small_progress) && evaluations < options.max_evaluations) {
            // Mesh vertices are nonsmooth. Confirm the stop using feasible projected coordinate probes
            // at several scales down to 1um, instead of treating a failed smooth model as convergence.
            Matrix poll_constraints(constraints.rows() + 2 * basis.rows(), dimensions);
            Vector poll_limits(poll_constraints.rows());
            poll_constraints.topRows(constraints.rows()) = constraints;
            poll_limits.head(constraints.rows()) = limits - constraints * variables;
            poll_constraints.middleRows(constraints.rows(), basis.rows()) = basis;
            poll_constraints.bottomRows(basis.rows()) = -basis;
            poll_limits.tail(2 * basis.rows()).setConstant(0.5);
            Vector best_variables = variables;
            double best_energy = current.energy;
            bool reliable = true;
            for (const Vector &poll_target : poll_targets) {
                for (const double sign : {-1., 1.}) {
                    if (evaluations >= options.max_evaluations) { reliable = false; break; }
                    Vector direction;
                    if (!quadratic_step(Matrix::Identity(dimensions, dimensions), sign * poll_target, poll_constraints, poll_limits, direction, &result.qp_failure)) {
                        reliable = false; break;
                    }
                    const double displacement = width * (basis * direction).lpNorm<Eigen::Infinity>();
                    if (displacement < 1e-10) continue;
                    double previous_fraction = -1.;
                    for (const double scale : {0.2, 0.1, 0.03, 0.01, 0.003, 0.001}) {
                        if (evaluations >= options.max_evaluations) { reliable = false; break; }
                        const double fraction = std::min(1., scale / displacement);
                        if (fraction == previous_fraction) continue;
                        previous_fraction = fraction;
                        const Vector candidate = variables + fraction * direction;
                        if (constraints.rows() > 0 && (constraints * candidate - limits).maxCoeff() > 1e-8) { reliable = false; break; }
                        const Evaluation probe = evaluate(contours, basis, candidate, width, options.length_weight, false);
                        ++evaluations;
                        if (!std::isfinite(probe.energy)) { reliable = false; break; }
                        // A microscopic cost change below the declared solver accuracy is not useful progress.
                        if (probe.energy < best_energy - 1e-8 * std::max(1., initial.energy)) {
                            best_variables = candidate;
                            best_energy = probe.energy;
                        }
                    }
                    if (!reliable) break;
                }
                if (!reliable) break;
            }
            if (reliable && best_energy < current.energy) {
                variables = best_variables;
                if (evaluations >= options.max_evaluations) { result.status = SeamChainStatus::Budget; return result; }
                current = evaluate(contours, basis, variables, width, options.length_weight, true);
                ++evaluations;
                result.trial_energy = current.energy;
                damping = 1e-3;
                accepted = true;
            } else if (reliable) {
                converged = true;
                result.resolution_converged = !stationary_model;
            }
        }
        if (converged) break;
        if (!accepted) {
            result.status = evaluations >= options.max_evaluations ? SeamChainStatus::Budget : SeamChainStatus::Failed;
            return result;
        }
    }
    if (!converged) { result.status = SeamChainStatus::Budget; return result; }
    // Compare against the exact original positions, including any tolerated reconstruction roundoff.
    if (current.energy >= initial.energy - 1e-8 * std::max(1., initial.energy) || current.length > result.initial_length + 1e-5) {
        result.status = SeamChainStatus::Unchanged;
        return result;
    }
    // Commit the whole chain only after convergence and geometry checks; failures never keep partial progress.
    for (size_t i = 0; i < points.size(); ++i)
        if (!current.positions[i].allFinite() || current.arcs[Eigen::Index(i)] < points[i].lower_arc - geometry_tolerance ||
            current.arcs[Eigen::Index(i)] > points[i].upper_arc + geometry_tolerance) return result;
    result.status = SeamChainStatus::Improved;
    result.positions = std::move(current.positions);
    result.arcs.assign(current.arcs.data(), current.arcs.data() + current.arcs.size());
    result.energy = current.energy;
    result.length = current.length;
    return result;
}
} // namespace Slic3r::SeamPlacerImpl
