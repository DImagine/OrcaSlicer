#include "SeamFeature.hpp"
#include "SeamChain.hpp"
#include "libslic3r/AABBTreeIndirect.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/libslic3r.h"
#include <admesh/stl.h>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace Slic3r::SeamPlacerImpl {
namespace {
constexpr double eps = 1e-5;
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

bool valid_options(const SeamFeatureOptions &o)
{
    const double positive[] = {o.corridor_widths, o.minimum_angle, o.minimum_edge_widths, o.offset_max_widths,
        o.geometry_tolerance, o.node_widths, o.horizontal_normal_min, o.capture_widths, o.transition_span,
        o.bridge_motion_widths, o.arc_chart_fraction, o.unbound_cost, o.solve_tolerance};
    for (double v : positive) if (!std::isfinite(v) || v <= 0.) return false;
    const double nonnegative[] = {o.trace_angle, o.offset_min_widths, o.offset_difference_widths,
        o.evidence_span_widths, o.route_margin, o.transition_bending, o.step_growth, o.step_allowance_widths,
        o.turn_allowance, o.support_allowance_widths};
    for (double v : nonnegative) if (!std::isfinite(v) || v < 0.) return false;
    return o.minimum_angle < 90. && o.trace_angle < 90. && o.arc_chart_fraction < .5 &&
        o.horizontal_normal_min < 1. && o.offset_min_widths < o.offset_max_widths && o.turn_allowance <= 180. &&
        o.evidence_layers >= 2 && o.evidence_window >= o.evidence_layers && o.evidence_window <= 101 &&
        o.max_iterations > 0 && o.max_chain_points > 1 && o.max_candidates > 0 && o.max_edge_queries > 0 &&
        o.bridge_repair_samples > 0 && o.bridge_repair_samples <= 256;
}

struct Loop {
    const SeamChainPoint *point = nullptr;
    std::vector<double> arcs;
    double initial = 0.;
    double total = 0.;
    double support_limit = 0.;
    bool prepare(const SeamChainPoint &p, const SeamFeatureOptions &o) {
        point = &p;
        if (p.vertices.size() < 3 || !p.initial.allFinite() || !std::isfinite(p.width) || p.width <= 0. ||
            !std::isfinite(p.anchor_weight) || p.anchor_weight < 0. || !std::isfinite(p.candidate_arc) ||
            std::isnan(p.lower_arc) || std::isnan(p.upper_arc) || p.lower_arc > p.upper_arc) return false;
        arcs.assign(1, 0.);
        double closest = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < p.vertices.size(); ++i) {
            const auto &a = p.vertices[i], &b = p.vertices[(i + 1) % p.vertices.size()];
            if (!a.position.allFinite() || !std::isfinite(a.overhang)) return false;
            const Vec2d d = b.position - a.position;
            const double length = d.norm();
            if (length > eps * .01) {
                const double t = std::clamp((p.initial.head<2>() - a.position).dot(d) / d.squaredNorm(), 0., 1.);
                const double distance = (p.initial.head<2>() - a.position - t * d).norm();
                if (distance < closest) { closest = distance; initial = arcs.back() + t * length; }
            }
            arcs.push_back(arcs.back() + length);
        }
        total = arcs.back();
        if (total <= eps || closest > eps) return false;
        initial += std::round((p.candidate_arc - initial) / total) * total;
        support_limit = support(p.candidate_arc) + o.support_allowance_widths * p.width;
        return initial >= p.lower_arc - eps && initial <= p.upper_arc + eps;
    }
    size_t edge(double &s) const {
        s = std::fmod(s, total); if (s < 0.) s += total;
        return std::min(point->vertices.size() - 1,
            size_t(std::upper_bound(arcs.begin(), arcs.end(), s) - arcs.begin() - 1));
    }
    Vec3d sample(double s) const {
        const size_t i = edge(s), j = (i + 1) % point->vertices.size();
        const double t = (s - arcs[i]) / (arcs[i + 1] - arcs[i]);
        const Vec2d p = (1. - t) * point->vertices[i].position + t * point->vertices[j].position;
        return Vec3d(p.x(), p.y(), point->initial.z());
    }
    double support(double s) const {
        const size_t i = edge(s), j = (i + 1) % point->vertices.size();
        const double t = (s - arcs[i]) / (arcs[i + 1] - arcs[i]);
        return (1. - t) * point->vertices[i].overhang + t * point->vertices[j].overhang;
    }
    double unwrap(double s) const { return s + std::round((initial - s) / total) * total; }
    std::pair<double, double> bounds(const SeamFeatureOptions &o) const {
        // Stay within one unambiguous arc chart, including on very small loops.
        const double radius = std::min(o.bridge_motion_widths * point->width, o.arc_chart_fraction * total);
        double low = std::max(point->lower_arc - initial, -radius);
        double high = std::min(point->upper_arc - initial, radius);
        const double limit = std::max(support_limit, support(initial));
        // Interpolated overhang is linear on each contour edge; restrict to the connected safe interval.
        for (const bool forward : {false, true}) {
            double cursor = initial, reach = 0.;
            const double wanted = forward ? high : -low;
            double current = support(cursor);
            for (size_t visited = 0; visited <= point->vertices.size() && reach < wanted; ++visited) {
                double wrapped = cursor + (forward ? eps * .001 : -eps * .001);
                const size_t i = edge(wrapped);
                const double boundary = forward ? arcs[i + 1] : arcs[i];
                double next = unwrap(boundary);
                if (forward && next <= cursor + eps * .001) next += total;
                if (!forward && next >= cursor - eps * .001) next -= total;
                const double step = std::min(wanted - reach, std::abs(next - cursor));
                const double end = cursor + (forward ? step : -step);
                const double other = support(end);
                if (other > limit + eps) {
                    reach += step * std::clamp((limit - current) / (other - current), 0., 1.);
                    break;
                }
                reach += step; cursor = end; current = other;
            }
            if (forward) high = reach; else low = -reach;
        }
        return {low, high};
    }
};

double point_segment_distance(const Vec3d &p, const Vec3d &a, const Vec3d &b)
{
    const Vec3d d = b - a;
    const double t = d.squaredNorm() > 0. ? std::clamp((p - a).dot(d) / d.squaredNorm(), 0., 1.) : 0.;
    return (p - a - t * d).norm();
}

double segment_distance(const Vec3d &a, const Vec3d &b, const Vec3d &p, const Vec3d &q)
{
    double result = std::min({point_segment_distance(a, p, q), point_segment_distance(b, p, q),
        point_segment_distance(p, a, b), point_segment_distance(q, a, b)});
    const Vec3d u = b - a, v = q - p, w = a - p;
    const double aa = u.squaredNorm(), bb = u.dot(v), cc = v.squaredNorm(), dd = u.dot(w), ee = v.dot(w);
    const double den = aa * cc - bb * bb;
    if (std::abs(den) > 1e-20) {
        const double s = (bb * ee - cc * dd) / den, t = (aa * ee - bb * dd) / den;
        if (s >= 0. && s <= 1. && t >= 0. && t <= 1.) result = std::min(result, (w + s * u - t * v).norm());
    }
    return result;
}

struct Graph {
    struct Edge { int a, b; std::vector<size_t> faces; double angle = -1.; bool core = false; };
    struct Box {
        size_t id; Eigen::AlignedBox<double, 3> bounds; Vec3d center;
        size_t idx() const { return id; }
        const Eigen::AlignedBox<double, 3> &bbox() const { return bounds; }
        const Vec3d &centroid() const { return center; }
    };
    const indexed_triangle_set &mesh;
    std::vector<Edge> edges;
    std::vector<Vec3d> normals;
    std::vector<std::vector<size_t>> incident;
    std::vector<bool> complete;
    std::vector<size_t> junction_degree;
    AABBTreeIndirect::Tree<3, double> tree;
    bool budget = false;
    Vec3d vertex(int i) const { return mesh.vertices[size_t(i)].cast<double>(); }
    double angle(size_t e) {
        auto &edge = edges[e];
        if (edge.angle < 0.) edge.angle = edge.faces.size() == 2 ?
            std::acos(std::clamp(normals[edge.faces[0]].dot(normals[edge.faces[1]]), -1., 1.)) * 180. / PI : 0.;
        return edge.angle;
    }
    Graph(const indexed_triangle_set &input, const std::vector<std::vector<SeamChainPoint>> &chains, const SeamFeatureOptions &o) : mesh(input) {
        normals.resize(mesh.indices.size(), Vec3d::Zero());
        std::map<std::pair<int, int>, std::vector<size_t>> adjacency;
        for (size_t i = 0; i < mesh.indices.size(); ++i) {
            const auto &f = mesh.indices[i];
            const Vec3d n = (vertex(f[1]) - vertex(f[0])).cross(vertex(f[2]) - vertex(f[0]));
            if (n.norm() <= 1e-12) continue;
            normals[i] = n.normalized();
            for (int j = 0; j < 3; ++j) adjacency[std::minmax(int(f[j]), int(f[(j + 1) % 3]))].push_back(i);
        }
        incident.resize(mesh.vertices.size()); complete.assign(mesh.vertices.size(), false);
        junction_degree.assign(mesh.vertices.size(), 0);
        std::vector<Box> boxes;
        for (const auto &[key, faces] : adjacency) {
            const size_t id = edges.size(); edges.push_back({key.first, key.second, faces});
            incident[size_t(key.first)].push_back(id); incident[size_t(key.second)].push_back(id);
            const Vec3d a = vertex(key.first), b = vertex(key.second);
            boxes.push_back({id, Eigen::AlignedBox<double, 3>(a.cwiseMin(b), a.cwiseMax(b)), .5 * (a + b)});
        }
        tree.build(std::move(boxes));
        size_t queries = 0;
        for (const auto &chain : chains) for (size_t i = 1; i < chain.size(); ++i) {
            const Vec3d p = chain[i - 1].initial, q = chain[i].initial;
            const double radius = o.corridor_widths * std::max(chain[i - 1].width, chain[i].width);
            const Vec3d extra = Vec3d::Constant(radius);
            const Eigen::AlignedBox<double, 3> box(p.cwiseMin(q) - extra, p.cwiseMax(q) + extra);
            AABBTreeIndirect::traverse(tree, AABBTreeIndirect::intersecting(box), [&](const auto &node) {
                const size_t e = node.idx; // The tree callback receives the leaf, not its primitive index.
                if (++queries > o.max_edge_queries) { budget = true; return false; }
                auto &edge = edges[e];
                if (edge.faces.size() == 2 && segment_distance(vertex(edge.a), vertex(edge.b), p, q) <= radius && angle(e) >= o.minimum_angle)
                    edge.core = true;
                return true;
            });
            if (budget) return; // Exhaustion never commits the partial graph.
        }
        for (auto &edge : edges) if (edge.core) complete[size_t(edge.a)] = complete[size_t(edge.b)] = true;
        for (size_t v = 0; v < complete.size(); ++v) if (complete[v])
            for (size_t e : incident[v]) if (angle(e) >= o.minimum_angle) ++junction_degree[v];
    }
    bool compatible(size_t a, size_t b, double low, double high, double width, const SeamFeatureOptions &o) const {
        if (a == size_t(-1) || b == size_t(-1) || a == b) return true;
        const auto &ea = edges[a], &eb = edges[b];
        for (int v : {ea.a, ea.b}) if ((v == eb.a || v == eb.b) && complete[size_t(v)] &&
            vertex(v).z() >= low - o.node_widths * width && vertex(v).z() <= high + o.node_widths * width) return true;
        return false;
    }
};

struct Choice { size_t edge, vertex; Vec3d target; double arc, distance, cost; };

std::vector<Choice> choices(Graph &graph, const Loop &loop, const SeamFeatureOptions &o)
{
    std::vector<Choice> result;
    const auto &p = *loop.point;
    for (size_t e = 0; e < graph.edges.size(); ++e) {
        const auto &edge = graph.edges[e];
        if (!edge.core || edge.faces.size() != 2) continue;
        const Vec3d a = graph.vertex(edge.a), b = graph.vertex(edge.b);
        if ((b - a).norm() < o.minimum_edge_widths * p.width || std::abs(b.z() - a.z()) < eps * .1) continue;
        const double node_radius = o.node_widths * p.width;
        if (p.initial.z() < std::min(a.z(), b.z()) - node_radius || p.initial.z() > std::max(a.z(), b.z()) + node_radius) continue;
        // Plane extensions are allowed only near a real junction, not at an isolated crease endpoint.
        if (p.initial.z() < std::min(a.z(), b.z()) || p.initial.z() > std::max(a.z(), b.z())) {
            const int endpoint = std::abs(p.initial.z() - a.z()) < std::abs(p.initial.z() - b.z()) ? edge.a : edge.b;
            if (graph.junction_degree[size_t(endpoint)] < 2) continue;
        }
        std::array<Vec2d, 2> normals, traces;
        std::array<double, 2> c;
        bool stable = true;
        for (size_t j = 0; j < 2; ++j) {
            const Vec3d &n = graph.normals[edge.faces[j]];
            const double h = n.head<2>().norm();
            if (h < o.horizontal_normal_min) { stable = false; break; }
            normals[j] = n.head<2>() / h; traces[j] = Vec2d(-normals[j].y(), normals[j].x());
            c[j] = (n.dot(a) - n.z() * p.initial.z()) / h;
        }
        if (!stable || std::abs(normals[0].x() * normals[1].y() - normals[0].y() * normals[1].x()) < std::sin(o.minimum_angle * PI / 180.)) continue;
        const double cosine = std::cos(o.trace_angle * PI / 180.);
        for (size_t v = 0; v < p.vertices.size(); ++v) {
            const auto &corner = p.vertices[v];
            const double distance = (corner.position - p.initial.head<2>()).norm();
            if (!corner.feature_allowed || distance > o.corridor_widths * p.width ||
                corner.overhang > loop.support_limit + eps * .01 || (p.fixed && distance > eps)) continue;
            size_t before = (v + p.vertices.size() - 1) % p.vertices.size(), after = (v + 1) % p.vertices.size();
            for (size_t n = 0; n < p.vertices.size() && (corner.position - p.vertices[before].position).norm() < eps * .01; ++n)
                before = (before + p.vertices.size() - 1) % p.vertices.size();
            for (size_t n = 0; n < p.vertices.size() && (corner.position - p.vertices[after].position).norm() < eps * .01; ++n)
                after = (after + 1) % p.vertices.size();
            const Vec2d in = (corner.position - p.vertices[before].position).normalized();
            const Vec2d out = (p.vertices[after].position - corner.position).normalized();
            if (!in.allFinite() || !out.allFinite() || in.dot(out) > std::cos(o.minimum_angle * PI / 180.)) continue;
            const bool match = (std::abs(in.dot(traces[0])) >= cosine && std::abs(out.dot(traces[1])) >= cosine) ||
                               (std::abs(in.dot(traces[1])) >= cosine && std::abs(out.dot(traces[0])) >= cosine);
            const double d0 = c[0] - normals[0].dot(corner.position), d1 = c[1] - normals[1].dot(corner.position);
            if (!match || std::min(d0, d1) < o.offset_min_widths * p.width || std::max(d0, d1) > o.offset_max_widths * p.width ||
                std::abs(d0 - d1) > o.offset_difference_widths * p.width + o.geometry_tolerance) continue;
            const double arc = loop.unwrap(loop.arcs[v]);
            if (arc < p.lower_arc - eps || arc > p.upper_arc + eps) continue;
            result.push_back({e, v, Vec3d(corner.position.x(), corner.position.y(), p.initial.z()), arc, distance,
                              std::pow(distance / p.width, 2.)});
        }
    }
    // Repeated sampling vertices may represent one geometric corner of the same carrier.
    std::sort(result.begin(), result.end(), [](const Choice &a, const Choice &b) { return std::tie(a.edge, a.arc) < std::tie(b.edge, b.arc); });
    result.erase(std::unique(result.begin(), result.end(), [](const Choice &a, const Choice &b) {
        return a.edge == b.edge && (a.target - b.target).norm() < eps;
    }), result.end());
    std::vector<Choice> unique;
    for (const Choice &choice : result) {
        // Two incident carriers can describe the same printed junction; disconnected copies cannot.
        const bool alias = std::any_of(unique.begin(), unique.end(), [&](const Choice &other) {
            return (choice.target - other.target).norm() < eps &&
                graph.compatible(choice.edge, other.edge, p.initial.z(), p.initial.z(), p.width, o);
        });
        if (!alias) unique.push_back(choice);
    }
    return unique;
}

bool box_solve(const Matrix &h, const Vector &low, const Vector &high, Vector &u, const SeamFeatureOptions &o)
{
    const size_t n = size_t(low.size());
    std::vector<int> state(n, 0);
    for (size_t i = 0; i < n; ++i) {
        if (low[Eigen::Index(i)] > high[Eigen::Index(i)]) return false;
        u[Eigen::Index(i)] = std::clamp(0., low[Eigen::Index(i)], high[Eigen::Index(i)]);
        if (high[Eigen::Index(i)] - low[Eigen::Index(i)] < o.solve_tolerance) state[i] = 2;
    }
    for (size_t iteration = 0; iteration < o.max_iterations; ++iteration) {
        std::vector<size_t> free;
        for (size_t i = 0; i < n; ++i) if (!state[i]) free.push_back(i);
        Vector proposed = u;
        if (!free.empty()) {
            Matrix sub(free.size(), free.size()); Vector rhs(free.size());
            for (size_t i = 0; i < free.size(); ++i) {
                rhs[Eigen::Index(i)] = 0.;
                for (size_t j = 0; j < n; ++j) if (state[j]) rhs[Eigen::Index(i)] -= h(Eigen::Index(free[i]), Eigen::Index(j)) * u[Eigen::Index(j)];
                for (size_t j = 0; j < free.size(); ++j) sub(Eigen::Index(i), Eigen::Index(j)) = h(Eigen::Index(free[i]), Eigen::Index(free[j]));
            }
            const Eigen::LDLT<Matrix> factor(sub);
            if (factor.info() != Eigen::Success || !factor.isPositive()) return false;
            const Vector solution = factor.solve(rhs);
            if (!solution.allFinite()) return false;
            for (size_t i = 0; i < free.size(); ++i) proposed[Eigen::Index(free[i])] = solution[Eigen::Index(i)];
        }
        double alpha = 1.; size_t hit = n; int side = 0;
        for (size_t i : free) {
            const double move = proposed[Eigen::Index(i)] - u[Eigen::Index(i)];
            if (move > 0. && proposed[Eigen::Index(i)] > high[Eigen::Index(i)]) {
                const double a = (high[Eigen::Index(i)] - u[Eigen::Index(i)]) / move;
                if (a < alpha) { alpha = a; hit = i; side = 1; }
            } else if (move < 0. && proposed[Eigen::Index(i)] < low[Eigen::Index(i)]) {
                const double a = (low[Eigen::Index(i)] - u[Eigen::Index(i)]) / move;
                if (a < alpha) { alpha = a; hit = i; side = -1; }
            }
        }
        u += std::max(0., alpha) * (proposed - u);
        if (hit < n) { state[hit] = side; u[Eigen::Index(hit)] = side > 0 ? high[Eigen::Index(hit)] : low[Eigen::Index(hit)]; continue; }
        const Vector gradient = h * u;
        size_t release = n; double violation = o.solve_tolerance;
        for (size_t i = 0; i < n; ++i) if (state[i] == -1 || state[i] == 1) {
            const double residual = state[i] == -1 ? -gradient[Eigen::Index(i)] : gradient[Eigen::Index(i)];
            if (residual > violation) { violation = residual; release = i; }
        }
        if (release == n) return true; // Feasible stationary solution of a strictly convex quadratic.
        state[release] = 0;
    }
    return false;
}

double step_limit(const SeamChainPoint &a, const SeamChainPoint &b, const SeamFeatureOptions &o, double carrier_step = 0.)
{
    // A smoothed first approximation may shorten the spacing of two true points on one inclined crease.
    const double old = std::max((a.initial - b.initial).norm(), carrier_step);
    return std::max(old * (1. + o.step_growth), old + o.step_allowance_widths * std::min(a.width, b.width));
}

double turn_angle(const Vec3d &a, const Vec3d &b, const Vec3d &c)
{
    const Vec3d u = b - a, v = c - b;
    return std::acos(std::clamp(u.dot(v) / (u.norm() * v.norm()), -1., 1.)) * 180. / PI;
}

using Interval = std::pair<double, double>;
std::vector<Interval> sphere_intervals(const Loop &loop, const Vec3d &center, double radius,
    const std::vector<Interval> &input)
{
    std::vector<Interval> output;
    const double r2 = radius * radius - std::pow(loop.point->initial.z() - center.z(), 2.);
    if (r2 < 0.) return output;
    for (const auto &[low, high] : input) for (size_t i = 0; i < loop.point->vertices.size(); ++i) {
        const double length = loop.arcs[i + 1] - loop.arcs[i];
        if (length <= eps * .01) continue;
        const Vec2d a = loop.point->vertices[i].position;
        const Vec2d direction = (loop.point->vertices[(i + 1) % loop.point->vertices.size()].position - a) / length;
        const Vec2d delta = a - center.head<2>();
        const double foot = -delta.dot(direction);
        const double margin = r2 - std::max(0., delta.squaredNorm() - foot * foot);
        if (margin < 0.) continue;
        const double half = std::sqrt(margin);
        const double first = std::max(0., foot - half), last = std::min(length, foot + half);
        if (first > last) continue;
        // The bounded chart covers less than a whole period; retain every intersecting contour edge.
        const double base = loop.arcs[i] - loop.initial;
        const int begin = int(std::floor((low - base - last) / loop.total));
        const int end = int(std::ceil((high - base - first) / loop.total));
        for (int cycle = begin; cycle <= end; ++cycle) {
            const double start = std::max(low, base + cycle * loop.total + first);
            const double finish = std::min(high, base + cycle * loop.total + last);
            if (start <= finish) output.emplace_back(start, finish);
        }
    }
    return output;
}

void repair_bridge(const Matrix &h, const Vector &low, const Vector &high, Vector &u,
    const std::vector<Loop> &loops, const std::vector<SeamChainPoint> &points,
    const std::vector<bool> &locked, const std::vector<bool> &geometric, const SeamFeatureOptions &o)
{
    std::vector<Vec3d> positions;
    for (size_t i = 0; i < points.size(); ++i) positions.push_back(loops[i].sample(loops[i].initial + u[Eigen::Index(i)]));
    for (size_t sweep = 0; sweep < o.bridge_repair_sweeps; ++sweep) {
        bool changed = false;
        for (size_t i = 0; i < points.size(); ++i) {
            if (locked[i]) continue;
            const bool bad_before = i && (positions[i] - positions[i - 1]).norm() > step_limit(points[i], points[i - 1], o) + eps;
            const bool bad_after = i + 1 < points.size() &&
                (positions[i] - positions[i + 1]).norm() > step_limit(points[i], points[i + 1], o) + eps;
            const size_t first_turn = i > 1 ? i - 1 : 1;
            const size_t last_turn = std::min(i + 1, points.size() - 2);
            bool bad_turn = false;
            if (o.bridge_repair_turns) for (size_t j = first_turn; j <= last_turn; ++j) {
                if (geometric[j - 1] && geometric[j] && geometric[j + 1]) continue;
                const double limit = turn_angle(points[j - 1].initial, points[j].initial, points[j + 1].initial) + o.turn_allowance;
                if (turn_angle(positions[j - 1], positions[j], positions[j + 1]) > limit + eps) bad_turn = true;
            }
            if (!bad_before && !bad_after && !bad_turn) continue;
            std::vector<Interval> allowed{{low[Eigen::Index(i)], high[Eigen::Index(i)]}};
            if (i) allowed = sphere_intervals(loops[i], positions[i - 1], step_limit(points[i], points[i - 1], o) + eps * .25, allowed);
            if (i + 1 < points.size()) allowed = sphere_intervals(loops[i], positions[i + 1], step_limit(points[i], points[i + 1], o) + eps * .25, allowed);
            if (allowed.size() > o.max_candidates) continue;
            // A one-coordinate quadratic minimum repairs both adjacent steps without moving any target.
            const double linear = h.row(Eigen::Index(i)).dot(u) - h(Eigen::Index(i), Eigen::Index(i)) * u[Eigen::Index(i)];
            const double optimum = -linear / h(Eigen::Index(i), Eigen::Index(i));
            double best = std::numeric_limits<double>::infinity(), chosen = u[Eigen::Index(i)];
            // Score every affected step and turn, so repairing one bend cannot hide another nearby violation.
            const auto violation = [&](double shift) {
                const Vec3d moved = loops[i].sample(loops[i].initial + shift);
                const auto at = [&](size_t j) -> Vec3d { return j == i ? moved : positions[j]; };
                double score = 0.;
                if (i) {
                    const double limit = step_limit(points[i], points[i - 1], o);
                    score = std::max(score, ((moved - positions[i - 1]).norm() - limit - eps) / std::max(limit, eps));
                }
                if (i + 1 < points.size()) {
                    const double limit = step_limit(points[i], points[i + 1], o);
                    score = std::max(score, ((moved - positions[i + 1]).norm() - limit - eps) / std::max(limit, eps));
                }
                for (size_t j = first_turn; j <= last_turn; ++j) {
                    if (geometric[j - 1] && geometric[j] && geometric[j + 1]) continue;
                    const double limit = turn_angle(points[j - 1].initial, points[j].initial, points[j + 1].initial) + o.turn_allowance;
                    score = std::max(score, (turn_angle(at(j - 1), at(j), at(j + 1)) - limit - eps) / 180.);
                }
                return score;
            };
            const double old_violation = o.bridge_repair_turns ? violation(chosen) : 0.;
            double best_violation = std::numeric_limits<double>::infinity();
            const auto consider = [&](double trial) {
                const double cost = .5 * h(Eigen::Index(i), Eigen::Index(i)) * trial * trial + linear * trial;
                const double score = o.bridge_repair_turns ? violation(trial) : 0.;
                if (score < best_violation - 1e-12 || (std::abs(score - best_violation) <= 1e-12 && cost < best)) {
                    best_violation = score; best = cost; chosen = trial;
                }
            };
            for (const auto &[start, finish] : allowed) {
                consider(std::clamp(optimum, start, finish));
                // Endpoints and bounded samples also find feasible turns missed by the quadratic optimum.
                if (o.bridge_repair_turns) for (size_t sample = 0; sample <= o.bridge_repair_samples; ++sample)
                    consider(start + (finish - start) * double(sample) / double(o.bridge_repair_samples));
            }
            if (std::isfinite(best) && (!o.bridge_repair_turns || best_violation <= old_violation + 1e-12) &&
                std::abs(chosen - u[Eigen::Index(i)]) > o.solve_tolerance) {
                u[Eigen::Index(i)] = chosen;
                positions[i] = loops[i].sample(loops[i].initial + chosen);
                changed = true;
            }
        }
        if (!changed) break; // Finite repair is followed by the complete step/turn validation, not presumed successful.
    }
}
} // namespace

bool seam_feature_continuous(const std::vector<SeamChainPoint> &original, const std::vector<Vec3d> &positions,
    const std::vector<bool> &geometric, const SeamFeatureOptions &o, size_t *failed,
    const std::vector<SeamFeatureBinding> *bindings)
{
    if (!valid_options(o) || positions.size() != original.size() || geometric.size() != original.size()) return false;
    if (bindings && bindings->size() != original.size()) return false;
    const auto fail = [&](size_t i) { if (failed) *failed = i; return false; };
    for (size_t i = 0; i < positions.size(); ++i) {
        if (!positions[i].allFinite() || std::abs(positions[i].z() - original[i].initial.z()) > eps ||
            (original[i].fixed && (positions[i] - original[i].initial).norm() > eps)) return fail(i);
        if (i == 0) continue;
        double carrier_step = 0.;
        if (o.carrier_reference_steps && bindings && geometric[i - 1] && geometric[i]) {
            const auto &a = (*bindings)[i - 1], &b = (*bindings)[i];
            if (a.edge != size_t(-1) && a.edge == b.edge && a.carrier_direction.allFinite() &&
                std::abs(a.carrier_direction.norm() - 1.) < eps && (a.carrier_direction - b.carrier_direction).norm() < eps &&
                std::abs(a.carrier_direction.z()) > 1e-12 && (a.target - positions[i - 1]).norm() <= eps &&
                (b.target - positions[i]).norm() <= eps)
                carrier_step = std::abs(positions[i].z() - positions[i - 1].z()) / std::abs(a.carrier_direction.z());
        }
        if ((positions[i] - positions[i - 1]).norm() > step_limit(original[i], original[i - 1], o, carrier_step) + eps) return fail(i);
        if (i >= 2 && !(geometric[i - 2] && geometric[i - 1] && geometric[i])) {
            const double old_turn = turn_angle(original[i - 2].initial, original[i - 1].initial, original[i].initial);
            if (turn_angle(positions[i - 2], positions[i - 1], positions[i]) > old_turn + o.turn_allowance + eps) return fail(i - 1);
        }
    }
    return true;
}

SeamFeatureResult accept_seam_feature_targets(const std::vector<SeamChainPoint> &points,
    std::vector<SeamFeatureBinding> bindings, const SeamFeatureOptions &o)
{
    SeamFeatureResult result;
    result.bindings = std::move(bindings);
    for (const auto &p : points) { result.positions.push_back(p.initial); result.locked.push_back(p.fixed); result.geometric.push_back(false); }
    if (!valid_options(o) || result.bindings.size() != points.size() || points.size() < 2 || points.size() > o.max_chain_points) return result;
    std::vector<Loop> loops(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        if (!loops[i].prepare(points[i], o) || (i && points[i].initial.z() <= points[i - 1].initial.z())) return result;
    const size_t n = points.size();
    Matrix h = Matrix::Zero(Eigen::Index(n), Eigen::Index(n));
    for (size_t i = 0; i < n; ++i) h(Eigen::Index(i), Eigen::Index(i)) = 1. + points[i].anchor_weight;
    const auto add = [&](const std::vector<size_t> &ids, const std::vector<double> &v, double weight) {
        for (size_t a = 0; a < ids.size(); ++a) for (size_t b = 0; b < ids.size(); ++b)
            h(Eigen::Index(ids[a]), Eigen::Index(ids[b])) += weight * v[a] * v[b];
    };
    for (size_t i = 1; i < n; ++i) {
        const double dz = points[i].initial.z() - points[i - 1].initial.z();
        add({i - 1, i}, {-1. / dz, 1. / dz}, o.transition_span * o.transition_span);
        if (i >= 2) {
            const double before = points[i - 1].initial.z() - points[i - 2].initial.z();
            add({i - 2, i - 1, i}, {1. / before, -1. / before - 1. / dz, 1. / dz},
                o.transition_bending * std::pow(o.transition_span, 4.) / std::pow(.5 * (dz + before), 2.));
        }
    }
    // Remove complete proposed runs on failure; never commit individual survivors of a broken bridge.
    for (size_t attempt = 0; attempt <= n; ++attempt) {
        Vector low(n), high(n), u = Vector::Zero(n);
        std::vector<bool> locked = result.locked;
        std::vector<bool> geometric(n, false);
        bool any = false;
        for (size_t i = 0; i < n; ++i) {
            const auto bounds = loops[i].bounds(o);
            low[Eigen::Index(i)] = bounds.first; high[Eigen::Index(i)] = bounds.second;
            if (points[i].fixed) low[Eigen::Index(i)] = high[Eigen::Index(i)] = 0.;
            auto &binding = result.bindings[i];
            if (binding.state == SeamFeatureState::Locked && points[i].fixed &&
                (binding.target - points[i].initial).norm() <= eps) geometric[i] = true;
            if (binding.state != SeamFeatureState::Proposed || points[i].fixed) continue;
            const double shift = binding.arc - loops[i].initial;
            if (binding.vertex >= points[i].vertices.size() || !points[i].vertices[binding.vertex].feature_allowed ||
                (loops[i].sample(binding.arc) - binding.target).norm() > eps || shift < bounds.first - eps || shift > bounds.second + eps ||
                loops[i].support(binding.arc) > loops[i].support_limit + eps) {
                binding.state = SeamFeatureState::TransitionRejected; continue;
            }
            low[Eigen::Index(i)] = high[Eigen::Index(i)] = shift; locked[i] = geometric[i] = true; any = true;
        }
        if (!any) return result;
        const bool solved = box_solve(h, low, high, u, o);
        if (solved) repair_bridge(h, low, high, u, loops, points, locked, geometric, o);
        std::vector<Vec3d> proposed(n);
        for (size_t i = 0; i < n; ++i) {
            proposed[i] = points[i].fixed ? points[i].initial : loops[i].sample(loops[i].initial + u[Eigen::Index(i)]);
            if (result.bindings[i].state == SeamFeatureState::Proposed && !points[i].fixed) proposed[i] = result.bindings[i].target;
        }
        size_t failed = 0;
        if (solved && seam_feature_continuous(points, proposed, geometric, o, &failed, &result.bindings)) {
            result.accepted = true; result.positions = std::move(proposed); result.locked = std::move(locked);
            result.geometric = std::move(geometric);
            for (size_t i = 0; i < n; ++i) if (result.bindings[i].state == SeamFeatureState::Proposed && !points[i].fixed)
                result.bindings[i].state = SeamFeatureState::Locked;
            return result;
        }
        if (!solved) {
            for (auto &binding : result.bindings) if (binding.state == SeamFeatureState::Proposed) binding.state = SeamFeatureState::Budget;
            return result;
        }
        size_t nearest = n; double distance = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < n; ++i) if (!points[i].fixed && result.bindings[i].state == SeamFeatureState::Proposed) {
            const double d = std::abs(points[i].initial.z() - points[failed].initial.z());
            if (d < distance) { nearest = i; distance = d; }
        }
        if (nearest == n) return result;
        size_t first = nearest, end = nearest + 1;
        while (first > 0 && !points[first - 1].fixed && result.bindings[first - 1].state == SeamFeatureState::Proposed) --first;
        while (end < n && !points[end].fixed && result.bindings[end].state == SeamFeatureState::Proposed) ++end;
        for (size_t i = first; i < end; ++i) result.bindings[i].state = SeamFeatureState::TransitionRejected;
    }
    return result;
}

std::vector<SeamFeatureResult> bind_seam_features(const indexed_triangle_set &mesh,
    const std::vector<std::vector<SeamChainPoint>> &chains, const SeamFeatureOptions &o)
{
    std::vector<SeamFeatureResult> results;
    if (!valid_options(o)) {
        for (const auto &chain : chains) results.push_back(accept_seam_feature_targets(chain, std::vector<SeamFeatureBinding>(chain.size()), o));
        return results;
    }
    Graph graph(mesh, chains, o);
    for (const auto &chain : chains) {
        std::vector<SeamFeatureBinding> bindings(chain.size());
        if (graph.budget || chain.size() > o.max_chain_points) {
            for (auto &b : bindings) b.state = SeamFeatureState::Budget;
            results.push_back(accept_seam_feature_targets(chain, std::move(bindings), o)); continue;
        }
        std::vector<std::vector<Choice>> states(chain.size());
        bool valid = chain.size() >= 2;
        for (size_t i = 0; i < chain.size() && valid; ++i) {
            Loop loop;
            valid = loop.prepare(chain[i], o);
            if (valid) states[i] = choices(graph, loop, o);
            if (states[i].size() > o.max_candidates) valid = false;
        }
        if (!valid) { results.push_back(accept_seam_feature_targets(chain, std::move(bindings), o)); continue; }
        for (auto &s : states) s.push_back({size_t(-1), size_t(-1), Vec3d::Zero(), 0., 0., o.unbound_cost});
        std::vector<std::vector<double>> forward(chain.size()), backward(chain.size());
        for (size_t i = 0; i < chain.size(); ++i) {
            forward[i].assign(states[i].size(), std::numeric_limits<double>::infinity());
            for (size_t b = 0; b < states[i].size(); ++b) {
                double best = i ? std::numeric_limits<double>::infinity() : 0.;
                if (i) for (size_t a = 0; a < states[i - 1].size(); ++a)
                    if (graph.compatible(states[i - 1][a].edge, states[i][b].edge, chain[i - 1].initial.z(), chain[i].initial.z(), chain[i].width, o))
                        best = std::min(best, forward[i - 1][a]);
                forward[i][b] = states[i][b].cost + best;
            }
        }
        backward.back().assign(states.back().size(), 0.);
        for (size_t i = chain.size() - 1; i-- > 0;) {
            backward[i].assign(states[i].size(), std::numeric_limits<double>::infinity());
            for (size_t a = 0; a < states[i].size(); ++a) for (size_t b = 0; b < states[i + 1].size(); ++b)
                if (graph.compatible(states[i][a].edge, states[i + 1][b].edge, chain[i].initial.z(), chain[i + 1].initial.z(), chain[i].width, o))
                    backward[i][a] = std::min(backward[i][a], states[i + 1][b].cost + backward[i + 1][b]);
        }
        const double best = *std::min_element(forward.back().begin(), forward.back().end());
        for (size_t i = 0; i < chain.size(); ++i) {
            std::vector<size_t> plausible, close;
            for (size_t j = 0; j < states[i].size(); ++j) {
                if (forward[i][j] + backward[i][j] <= best + o.route_margin) plausible.push_back(j);
                if (states[i][j].edge != size_t(-1) && states[i][j].distance <= o.capture_widths * chain[i].width) close.push_back(j);
            }
            if (close.size() > 1 || plausible.size() > 1) { bindings[i].state = SeamFeatureState::Ambiguous; continue; }
            if (close.size() != 1 || plausible.size() != 1 || close[0] != plausible[0]) continue;
            const Choice &choice = states[i][close[0]];
            size_t evidence = 0; double low = chain[i].initial.z(), high = low;
            const size_t half = o.evidence_window / 2;
            for (size_t j = i > half ? i - half : 0; j < std::min(chain.size(), i + half + 1); ++j)
                if (std::any_of(states[j].begin(), states[j].end(), [&](const Choice &s) { return s.edge == choice.edge; })) {
                    ++evidence; low = std::min(low, chain[j].initial.z()); high = std::max(high, chain[j].initial.z());
                }
            if (evidence < o.evidence_layers || high - low + o.solve_tolerance < o.evidence_span_widths * chain[i].width) {
                bindings[i].state = SeamFeatureState::InsufficientEvidence; continue;
            }
            bindings[i] = {chain[i].fixed ? SeamFeatureState::Locked : SeamFeatureState::Proposed,
                           choice.edge, choice.vertex, choice.target, choice.arc};
            const auto &carrier = graph.edges[choice.edge];
            bindings[i].carrier_direction = (graph.vertex(carrier.b) - graph.vertex(carrier.a)).normalized();
        }
        results.push_back(accept_seam_feature_targets(chain, std::move(bindings), o));
    }
    return results;
}
} // namespace Slic3r::SeamPlacerImpl
