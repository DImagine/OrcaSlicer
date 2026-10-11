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

double seam_chain_anchor_weight(double base_weight, double legacy_factor, double angle_degrees,
                               bool painted, const SeamChainAngleOptions &options)
{
    // Legacy mode reproduces the adapter's old arithmetic, including its different painted factor.
    if (options.curve == SeamChainAngleCurve::Legacy)
        return painted ? base_weight + options.painted_center_weight : base_weight * legacy_factor;
    const double x = std::clamp((std::abs(angle_degrees) - options.minimum_degrees) /
                               (options.full_degrees - options.minimum_degrees), 0., 1.);
    const double response = options.curve == SeamChainAngleCurve::Linear ? x :
                            options.curve == SeamChainAngleCurve::Cubic ? x * x * x : x * x * (3. - 2. * x);
    // Full soft weight does not fix a point: the original 55-degree mask is computed separately.
    return base_weight * response + (painted ? options.painted_center_weight : 0.);
}

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
    double preference_cost = 0.;
    double length_cost = 0.;
    double projection_cost = 0.;
    Vector gradient;
    Matrix hessian;
    std::vector<Vec3d> positions;
    Vector arcs;
};

Evaluation evaluate(const std::vector<Contour> &contours, const Matrix &basis, const Vector &variables,
                    double width, const SeamChainOptions &options, bool derivatives)
{
    const double lambda = options.length_weight, vertical_weight = options.vertical_weight;
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
        const double preference = height_weight * (anchor_factor * distance * distance + q_residual * q_residual + o_residual * o_residual);
        value.preference_cost += preference;
        value.energy += preference;
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
            const double length_cost = lambda / width * delta.head<2>().squaredNorm() / (length + delta.z());
            value.length_cost += length_cost;
            value.energy += length_cost;
            if (derivatives) {
                const Matrix difference = jacobian - last_jacobian;
                value.gradient += lambda / (width * length) * difference.transpose() * delta;
                // A positive quadratic majorizer of length avoids an indefinite Newton system.
                value.hessian += lambda / (width * length) * difference.transpose() * difference;
            }
            const SeamChainPoint &previous = *contours[i - 1].point;
            const auto free_point = [](const SeamChainPoint &p) {
                return !p.fixed && !std::isfinite(p.lower_arc) && !std::isfinite(p.upper_arc);
            };
            if (options.direction_layers == 0 && vertical_weight > 0. && free_point(previous) && free_point(point) &&
                previous.surface_confidence > 0. && point.surface_confidence > 0. &&
                previous.surface_normal.dot(point.surface_normal) > std::cos(options.sharp_surface_angle * PI / 180.)) {
                // Project the seam and Z onto the same frozen tangent plane: their sideways axis is n cross Z.
                const Vec3d normal = (previous.surface_normal + point.surface_normal).normalized();
                Vec2d side(normal.y(), -normal.x());
                const double horizontal = side.squaredNorm();
                if (horizontal > 1e-8) {
                    side.normalize();
                    const double confidence = std::min(previous.surface_confidence, point.surface_confidence);
                    // Near-horizontal surfaces have almost no projected Z and cannot define reliable verticality.
                    const double factor = vertical_weight * confidence * horizontal / (width * delta.z());
                    const double residual = side.dot(delta.head<2>());
                    const double projection = factor * residual * residual;
                    value.projection_cost += projection;
                    value.energy += projection;
                    if (derivatives) {
                        const Vector slope = (jacobian - last_jacobian).topRows(2).transpose() * side;
                        value.gradient += 2. * factor * residual * slope;
                        value.hessian += 2. * factor * slope * slope.transpose();
                    }
                }
            }
        }
        last_jacobian = jacobian;
    }
    if (options.direction_layers == 3 && vertical_weight > 0.) for (size_t i = 0; i < count; ++i) {
        const auto &point = *contours[i].point;
        if (point.fixed || point.surface_confidence <= 0.) continue;
        std::vector<size_t> members;
        double mean_z = 0.;
        for (size_t j = i > 0 ? i - 1 : i; j <= std::min(i + 1, count - 1); ++j) {
            const auto &neighbor = *contours[j].point;
            // Unreliable normals remove samples from this cost, never their contour motion variables.
            if (neighbor.surface_confidence <= 0. || neighbor.surface_normal.dot(point.surface_normal) <
                std::cos(options.sharp_surface_angle * PI / 180.)) continue;
            members.push_back(j);
            mean_z += neighbor.initial.z();
        }
        if (members.size() < 2) continue;
        mean_z /= double(members.size());
        double variance = 0.;
        for (size_t j : members) variance += std::pow(contours[j].point->initial.z() - mean_z, 2.);
        if (variance < 1e-14) continue;
        Vec2d side(point.surface_normal.y(), -point.surface_normal.x());
        const double horizontal = side.squaredNorm();
        if (horizontal < 1e-8) continue;
        side.normalize();
        double residual = 0.;
        Vector slope = Vector::Zero(dimensions);
        for (size_t j : members) {
            const double coefficient = (contours[j].point->initial.z() - mean_z) / variance;
            residual += coefficient * side.dot(value.positions[j].head<2>());
            if (derivatives) slope += coefficient * width *
                side.dot(sample(contours[j], value.arcs[Eigen::Index(j)]).tangent.head<2>()) * basis.row(Eigen::Index(j)).transpose();
        }
        double height = 0.;
        if (i > 0) height += 0.5 * (point.initial.z() - contours[i - 1].point->initial.z());
        if (i + 1 < count) height += 0.5 * (contours[i + 1].point->initial.z() - point.initial.z());
        const double factor = vertical_weight * point.surface_confidence * horizontal * height / width;
        value.projection_cost += factor * residual * residual;
        value.energy += factor * residual * residual;
        if (derivatives) {
            value.gradient += 2. * factor * residual * slope;
            value.hessian += 2. * factor * slope * slope.transpose();
        }
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

void prepare_seam_projection(std::vector<SeamChainPoint> &points, const SeamChainOptions &options)
{
    if (points.empty()) return;
    std::vector<double> arc(points.size(), 0.);
    for (size_t i = 1; i < points.size(); ++i)
        arc[i] = arc[i - 1] + (points[i].initial - points[i - 1].initial).norm();
    const bool anchored = std::any_of(points.begin(), points.end(), [](const SeamChainPoint &p) {
        return p.fixed || std::isfinite(p.lower_arc) || std::isfinite(p.upper_arc);
    });
    const auto valid_frame = [&](size_t i) {
        const auto &p = points[i];
        return p.surface_confidence > 0. && p.surface_normal.allFinite() &&
            std::abs(p.surface_normal.norm() - 1.) < 1e-5 &&
            std::abs(p.surface_normal.z()) <= std::sin(options.max_face_tilt * PI / 180.);
    };
    const double agreement = std::cos(options.sharp_surface_angle * PI / 180.);
    std::vector<bool> held(points.size(), false);
    if (options.projection_gate) for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].fixed) continue;
        if (!valid_frame(i)) { points[i].surface_confidence = 0.; continue; }
        if (options.direction_layers == 3) {
            Vec3d trend = Vec3d::Zero();
            std::vector<size_t> members;
            double mean_z = 0.;
            for (size_t j = i > 0 ? i - 1 : i; j <= std::min(i + 1, points.size() - 1); ++j) {
                if (!valid_frame(j) || points[j].surface_normal.dot(points[i].surface_normal) < agreement) continue;
                members.push_back(j);
                mean_z += points[j].initial.z();
            }
            if (members.size() < 2) { points[i].surface_confidence = 0.; continue; }
            mean_z /= double(members.size());
            for (size_t j : members) trend += (points[j].initial.z() - mean_z) * (points[j].initial - points[i].initial);
            const Vec3d &normal = points[i].surface_normal;
            const Vec3d side = Vec3d(-normal.y(), normal.x(), 0.).normalized();
            const Vec3d vertical = (Vec3d::UnitZ() - normal.z() * normal).normalized();
            if (trend.norm() > 1e-12)
                points[i].projection_angle = std::atan2(std::abs(trend.dot(side)), std::abs(trend.dot(vertical))) * 180. / PI;
            held[i] = points[i].projection_angle > options.max_projection_angle;
            continue;
        }
        size_t first = i, last = i;
        // Clip the analysis window at an unreliable surface or an existing sharp anchor.
        while (first > 0 && arc[i] - arc[first] < 0.5 * options.projection_span &&
               valid_frame(first - 1) && !points[first - 1].fixed &&
               points[first].surface_normal.dot(points[first - 1].surface_normal) >= agreement) --first;
        while (last + 1 < points.size() && arc[last] - arc[i] < 0.5 * options.projection_span &&
               valid_frame(last + 1) && !points[last + 1].fixed &&
               points[last].surface_normal.dot(points[last + 1].surface_normal) >= agreement) ++last;
        const double low = std::max(arc[first], arc[i] - 0.5 * options.projection_span);
        const double high = std::min(arc[last], arc[i] + 0.5 * options.projection_span);
        const double center = 0.5 * (low + high);
        Vec3d trend = Vec3d::Zero();
        // Integrate the linear-regression moment exactly on the clipped seam polyline.
        for (size_t j = first; j < last; ++j) {
            const double a = std::max(low, arc[j]), b = std::min(high, arc[j + 1]);
            if (b <= a || arc[j + 1] <= arc[j]) continue;
            const Vec3d slope = (points[j + 1].initial - points[j].initial) / (arc[j + 1] - arc[j]);
            const double middle = 0.5 * (a + b), length = b - a;
            const Vec3d position = points[j].initial + (middle - arc[j]) * slope;
            trend += length * (middle - center) * (position - points[i].initial) +
                     (length * length * length / 12.) * slope;
        }
        const Vec3d &normal = points[i].surface_normal;
        const Vec3d side = Vec3d(-normal.y(), normal.x(), 0.).normalized();
        const Vec3d vertical = (Vec3d::UnitZ() - normal.z() * normal).normalized();
        if (trend.norm() > 1e-12)
            points[i].projection_angle = std::atan2(std::abs(trend.dot(side)), std::abs(trend.dot(vertical))) * 180. / PI;
        held[i] = points[i].projection_angle > options.max_projection_angle;
    }
    if (options.hold_free_center && !anchored) {
        // Pick an existing layer point nearest half the 3D arc length, before eligibility freezes invent anchors.
        const auto middle = std::min_element(arc.begin(), arc.end(), [&](double a, double b) {
            return std::abs(a - 0.5 * arc.back()) < std::abs(b - 0.5 * arc.back());
        });
        points[size_t(middle - arc.begin())].fixed = true;
    }
    for (size_t i = 0; i < points.size(); ++i) if (held[i] && !points[i].fixed) {
        points[i].projection_hold = points[i].fixed = true;
    }
}

SeamChainResult optimize_seam_chain(const std::vector<SeamChainPoint> &points, const SeamChainOptions &options)
{
    SeamChainResult result;
    for (size_t i = 0; i < points.size(); ++i) {
        result.positions.push_back(points[i].initial);
        if (i > 0) result.initial_length += (points[i].initial - points[i - 1].initial).norm();
    }
    result.length = result.initial_length;
    if (options.length_weight == 0.) { result.status = SeamChainStatus::Disabled; return result; }
    if (points.size() < 2 || options.segments == 0 || !std::isfinite(options.length_weight) || options.length_weight < 0. ||
        !std::isfinite(options.vertical_weight) || options.vertical_weight < 0. ||
        !std::isfinite(options.projection_span) || options.projection_span <= 0. ||
        !std::isfinite(options.relative_length_growth) || options.relative_length_growth < 0. ||
        !std::isfinite(options.max_motion_widths) || options.max_motion_widths <= 0. ||
        !std::isfinite(options.boundary_taper_fraction) || options.boundary_taper_fraction <= 0. || options.boundary_taper_fraction > 0.5 ||
        !std::isfinite(options.max_projection_angle) || options.max_projection_angle < 0. || options.max_projection_angle > 90. ||
        !std::isfinite(options.max_face_tilt) || options.max_face_tilt < 0. || options.max_face_tilt >= 90. ||
        !std::isfinite(options.sharp_surface_angle) || options.sharp_surface_angle < 0. || options.sharp_surface_angle > 180. ||
        (options.direction_layers != 0 && options.direction_layers != 3)) return result;
    const bool projected_pass = options.vertical_weight > 0. || options.geometric_probe;
    const double length_limit = result.initial_length * (1. + (projected_pass ? options.relative_length_growth : 0.)) + geometry_tolerance;
    const auto motion_radius = [&](size_t i) {
        return std::min(options.max_motion_widths * points[i].width, points[i].surface_radius);
    };
    std::vector<Contour> contours(points.size());
    std::vector<double> widths;
    std::vector<size_t> fixed;
    for (size_t i = 0; i < points.size(); ++i) {
        if (!prepare(points[i], contours[i]) || (i > 0 && points[i].initial.z() <= points[i - 1].initial.z())) return result;
        if (!points[i].surface_normal.allFinite() || !std::isfinite(points[i].surface_confidence) ||
            points[i].surface_confidence < 0. || points[i].surface_confidence > 1.) return result;
        if (points[i].surface_confidence > 0. && std::abs(points[i].surface_normal.norm() - 1.) > 1e-5) return result;
        if (std::isnan(points[i].surface_radius) || points[i].surface_radius < 0.) return result;
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
    const bool patches = std::any_of(points.begin(), points.end(), [](const SeamChainPoint &p) {
        return p.surface_confidence > 0.;
    });
    if (options.projection_gate) {
        // Independent smooth deformation fields retain freedom when ineligible rows would kill a global spline.
        basis.resize(Eigen::Index(points.size()), 0);
        for (size_t first = 0; first < points.size();) {
            if (points[first].fixed) { ++first; continue; }
            size_t end = first + 1;
            while (end < points.size() && !points[end].fixed) ++end;
            const double low = points[first > 0 ? first - 1 : first].initial.z();
            const double high = points[end < points.size() ? end : end - 1].initial.z();
            const double transition = std::min(0.5 * options.projection_span, options.boundary_taper_fraction * (high - low));
            const size_t pieces = std::max(size_t(1), size_t(std::ceil((high - low) / options.projection_span)));
            const Eigen::Index column = basis.cols();
            basis.conservativeResize(Eigen::NoChange, column + Eigen::Index(pieces + 1));
            basis.rightCols(Eigen::Index(pieces + 1)).setZero();
            for (size_t i = first; i < end; ++i) {
                const double z = high > low ? (points[i].initial.z() - low) * double(pieces) / (high - low) : 0.;
                double taper = 1.;
                const auto smooth = [](double x) { x = std::clamp(x, 0., 1.); return x * x * (3. - 2. * x); };
                if (first > 0) taper *= smooth((points[i].initial.z() - low) / transition);
                if (end < points.size()) taper *= smooth((high - points[i].initial.z()) / transition);
                for (int j = int(std::floor(z)) - 1; j < int(std::floor(z)) + 3; ++j)
                    basis(Eigen::Index(i), column + std::clamp(j, 0, int(pieces))) +=
                        taper * Geometry::CubicBSplineKernel<double>::kernel(double(j) - z);
            }
            first = end;
        }
    }
    if (options.geometric_probe && patches) {
        const auto free_point = [&](size_t i) {
            const auto &p = points[i];
            return !p.fixed && !std::isfinite(p.lower_arc) && !std::isfinite(p.upper_arc) && p.surface_confidence > 0.;
        };
        for (size_t first = 0; first < points.size();) {
            if (!free_point(first)) { ++first; continue; }
            size_t end = first + 1;
            while (end < points.size() && free_point(end)) ++end;
            Vector profile = Vector::Zero(basis.rows());
            const double low = points[first > 0 ? first - 1 : first].initial.z();
            const double high = points[end < points.size() ? end : end - 1].initial.z();
            const double transition = std::min(0.5 * options.projection_span, options.boundary_taper_fraction * (high - low));
            const auto anchor = [&](size_t i) {
                return (points[i].fixed && !points[i].projection_hold) || std::isfinite(points[i].lower_arc) || std::isfinite(points[i].upper_arc);
            };
            // Loss of normal confidence only clips a patch; it must not invent a positional anchor.
            const bool left = first > 0 && anchor(first - 1);
            const bool right = end < points.size() && anchor(end);
            for (size_t i = first; i < end; ++i) {
                if (points[i].surface_confidence <= 0.) continue;
                // This geometric seed is a probe, never a fixed endpoint or a positional penalty.
                Vec2d target = points[(first + end - 1) / 2].initial.head<2>();
                if (left && right) {
                    const double fraction = (points[i].initial.z() - points[first - 1].initial.z()) /
                        (points[end].initial.z() - points[first - 1].initial.z());
                    target = (1. - fraction) * points[first - 1].initial.head<2>() + fraction * points[end].initial.head<2>();
                } else if (left) target = points[first - 1].initial.head<2>();
                else if (right) target = points[end].initial.head<2>();
                // Aim sideways in this fixed plane; normal-depth changes must not attract the probe to another wall.
                Vec2d side(points[i].surface_normal.y(), -points[i].surface_normal.x());
                if (side.norm() < 1e-6) continue;
                side.normalize();
                const Vec2d initial_xy = points[i].initial.head<2>();
                target = initial_xy + side.dot(target - initial_xy) * side;
                double nearest = std::numeric_limits<double>::infinity();
                double target_arc = contours[i].initial_arc;
                for (size_t j = 0; j < points[i].vertices.size(); ++j) {
                    const Vec2d a = points[i].vertices[j].position;
                    const Vec2d edge = points[i].vertices[(j + 1) % points[i].vertices.size()].position - a;
                    if (edge.squaredNorm() == 0.) continue;
                    const double t = std::clamp((target - a).dot(edge) / edge.squaredNorm(), 0., 1.);
                    const double distance = (target - a - t * edge).squaredNorm();
                    if (distance < nearest) { nearest = distance; target_arc = contours[i].arc[j] + t * edge.norm(); }
                }
                const double shift = std::remainder(target_arc - contours[i].initial_arc, contours[i].total);
                // Stay on the same local surface chart rather than crossing to an opposite thin wall.
                if (std::abs(shift) <= motion_radius(i)) {
                    double taper = 1.;
                    const auto smooth = [](double x) { x = std::clamp(x, 0., 1.); return x * x * (3. - 2. * x); };
                    // Eligibility boundaries retain a smooth displacement transition, rather than a one-layer jump.
                    if (options.projection_gate && first > 0 && points[first - 1].fixed)
                        taper *= smooth((points[i].initial.z() - points[first - 1].initial.z()) / transition);
                    if (options.projection_gate && end < points.size() && points[end].fixed)
                        taper *= smooth((points[end].initial.z() - points[i].initial.z()) / transition);
                    profile[Eigen::Index(i)] = taper * shift;
                }
            }
            const auto extend = [&](size_t neighbor, size_t reliable) {
                const auto &p = points[neighbor];
                if (!p.fixed && p.surface_confidence == 0. && !std::isfinite(p.lower_arc) && !std::isfinite(p.upper_arc))
                    profile[Eigen::Index(neighbor)] = std::clamp(profile[Eigen::Index(reliable)], -motion_radius(neighbor), motion_radius(neighbor));
            };
            // A missing normal separates independently searchable routes, but its point can follow either neighbor.
            if (first > 0) extend(first - 1, first);
            if (end < points.size()) extend(end, end - 1);
            first = end;
            const double amplitude = profile.lpNorm<Eigen::Infinity>();
            if (amplitude > 1e-8) {
                // Give each free continuation its own probe; unrelated runs must not share a forced amplitude.
                profile /= amplitude;
                const Vector residual = profile - basis * basis.colPivHouseholderQr().solve(profile);
                if (residual.norm() > 1e-6 * profile.norm()) {
                    basis.conservativeResize(Eigen::NoChange, basis.cols() + 1);
                    basis.col(basis.cols() - 1) = profile;
                }
            }
        }
    }
    if (options.projection_gate && basis.cols() > 0) {
        // Very short eligible runs may produce dependent columns; remove them before QP factorization.
        Eigen::ColPivHouseholderQR<Matrix> factor(basis);
        basis = (basis * factor.colsPermutation()).leftCols(factor.rank()).eval();
    }
    if (!fixed.empty()) {
        if (basis.cols() == 0) { result.status = SeamChainStatus::NoFreedom; return result; }
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
        if (patches && (options.vertical_weight > 0. || options.geometric_probe)) {
            // Total second-pass motion stays within the fixed patch; this is separate from painted bounds.
            rows.push_back(basis.row(Eigen::Index(i)).transpose());
            upper.push_back(motion_radius(i) / width);
            rows.push_back(-basis.row(Eigen::Index(i)).transpose());
            upper.push_back(motion_radius(i) / width);
        }
    }
    Matrix constraints(Eigen::Index(rows.size()), dimensions);
    Vector limits(Eigen::Index(rows.size()));
    for (size_t i = 0; i < rows.size(); ++i) { constraints.row(Eigen::Index(i)) = rows[i].transpose(); limits[Eigen::Index(i)] = upper[i]; }
    Vector variables = Vector::Zero(dimensions);
    if (options.max_evaluations == 0) { result.status = SeamChainStatus::Budget; return result; }
    Evaluation initial = evaluate(contours, basis, variables, width, options, true);
    Evaluation current = initial;
    result.initial_energy = result.energy = initial.energy;
    result.preference_cost = initial.preference_cost;
    result.length_cost = initial.length_cost;
    result.projection_cost = initial.projection_cost;
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
                Evaluation trial = evaluate(contours, basis, trial_variables, width, options, false);
                if (projected_pass && trial.length > length_limit) { ++result.length_rejections; continue; }
                if (!std::isfinite(trial.energy) || trial.energy >= current.energy - 1e-14 * std::max(1., current.energy)) continue;
                small_progress = current.energy - trial.energy <= 1e-8 * std::max(1., initial.energy) &&
                                 std::ldexp(movement, -int(backtrack)) <= 0.001;
                variables = trial_variables;
                if (evaluations >= options.max_evaluations) { result.status = SeamChainStatus::Budget; return result; }
                current = evaluate(contours, basis, variables, width, options, true);
                result.trial_energy = current.energy;
                result.trial_length = current.length;
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
                        const Evaluation probe = evaluate(contours, basis, candidate, width, options, false);
                        ++evaluations;
                        if (!std::isfinite(probe.energy)) { reliable = false; break; }
                        if (projected_pass && probe.length > length_limit) { ++result.length_rejections; continue; }
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
                current = evaluate(contours, basis, variables, width, options, true);
                ++evaluations;
                result.trial_energy = current.energy;
                result.trial_length = current.length;
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
    if (current.energy >= initial.energy - 1e-8 * std::max(1., initial.energy) || current.length > length_limit) {
        result.status = SeamChainStatus::Unchanged;
        return result;
    }
    // Commit the whole chain only after convergence and geometry checks; failures never keep partial progress.
    for (size_t i = 0; i < points.size(); ++i)
        if (!current.positions[i].allFinite() || current.arcs[Eigen::Index(i)] < points[i].lower_arc - geometry_tolerance ||
            current.arcs[Eigen::Index(i)] > points[i].upper_arc + geometry_tolerance ||
            (patches && (options.vertical_weight > 0. || options.geometric_probe) &&
             std::abs(current.arcs[Eigen::Index(i)] - contours[i].initial_arc) >
               motion_radius(i) + geometry_tolerance)) return result;
    result.status = SeamChainStatus::Improved;
    result.positions = std::move(current.positions);
    result.arcs.assign(current.arcs.data(), current.arcs.data() + current.arcs.size());
    result.energy = current.energy;
    result.preference_cost = current.preference_cost;
    result.length_cost = current.length_cost;
    result.projection_cost = current.projection_cost;
    result.length = current.length;
    return result;
}
} // namespace Slic3r::SeamPlacerImpl
