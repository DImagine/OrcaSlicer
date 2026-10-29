#include "PreciseSeamRulesGeometry.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>
#include "PreciseSeamInternal.hpp"
#include "libslic3r/AABBTreeLines.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Line.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"

// Precise Seam with layer rules: the seam point of a perimeter for the modifiers of an object that uses rules,
// with their effective types, layers and magnets. Objects without rules keep the drivers of PreciseSeam.cpp.
// Design: docs/HLSD/precise-seam.md

namespace Slic3r {
namespace PreciseSeam {

using SeamPlacerImpl::EnforcedBlockedSeamPoint;

// --- Magnet ----------------------------------------------------------------------------------------------------

// A rear point goes first, then a left one: the last tie-break of strong segments and magnet points.
// Slice coordinates already include object rotation and retain bed axes; centering and instance
// translation cannot change this order.
static bool rear_then_left(const Point &a, const Point &b)
{
    return a.y() > b.y() || (a.y() == b.y() && a.x() < b.x());
}

// The nearest point of segment [a, b] to `point` and its parameter on the segment.
static std::pair<Vec2d, double> project_on_segment(const Vec2d &point, const Vec2d &a, const Vec2d &b)
{
    const Vec2d  ab     = b - a;
    const double length = ab.squaredNorm();
    const double t      = length > 0. ? std::clamp((point - a).dot(ab) / length, 0., 1.) : 0.;
    return {a + t * ab, t};
}

// Positions along a perimeter as arc lengths from vertex 0.
struct ArcLengths {
    std::vector<double> start; // Of each vertex, and the total length last.

    explicit ArcLengths(const Lines &lines) : start(lines.size() + 1, 0.)
    {
        for (size_t edge = 0; edge < lines.size(); ++edge)
            start[edge + 1] = start[edge] + lines[edge].length();
    }
    double total() const { return start.back(); }
    double along(const PerimeterPosition &p) const
    {
        return start[p.edge_index] + p.parameter * (start[p.edge_index + 1] - start[p.edge_index]);
    }
    // The length forward from one position to another, around the start.
    double forward(double from, double to) const { return to >= from ? to - from : to - from + total(); }
    // The position at an arc length below the total.
    PerimeterPosition position(double at) const
    {
        const size_t edge   = std::min<size_t>(std::upper_bound(start.begin(), start.end(), at) - start.begin(), start.size() - 1) - 1;
        const double length = start[edge + 1] - start[edge];
        return {edge, length > 0. ? std::clamp((at - start[edge]) / length, 0., 1.) : 0.};
    }
    // The point at an arc length.
    Vec2d point(const Lines &lines, double at) const
    {
        const size_t edge = std::min<size_t>(std::upper_bound(start.begin(), start.end(), at) - start.begin(), lines.size()) - 1;
        const double length = start[edge + 1] - start[edge];
        const double t      = length > 0. ? std::clamp((at - start[edge]) / length, 0., 1.) : 0.;
        return lines[edge].a.cast<double>() + t * (lines[edge].b - lines[edge].a).cast<double>();
    }
};

// The tip of a corner: the vertex strictly inside the part `span` long forward from `from` that is farthest
// from `centre`, when it lies deeper than `ends` (the depth of the part's ends) by more than the tolerance.
// The distance along a segment is convex, so no point between vertices lies deeper. Equally deep vertices go
// to the one nearest the perimeter centroid, then to the rear and left one.
static std::optional<size_t> deepest_vertex(const Polygon &perimeter, const ArcLengths &arcs, double from, double span,
                                            const Vec2d &centre, double ends, double tolerance)
{
    std::vector<std::pair<size_t, double>> inside;
    double                                 deepest = 0.;
    for (size_t k = 0; k < perimeter.size(); ++k) {
        const double offset = arcs.forward(from, arcs.start[k]);
        if (offset <= 0. || offset >= span)
            continue;
        const double distance = (perimeter[k].cast<double>() - centre).norm();
        inside.emplace_back(k, distance);
        deepest = std::max(deepest, distance);
    }
    if (inside.empty() || deepest <= ends + tolerance)
        return std::nullopt;
    std::vector<size_t> candidates;
    for (const auto &[k, distance] : inside)
        if (distance >= deepest - tolerance)
            candidates.push_back(k);
    if (perimeter.area() != 0.) {
        const Vec2d centroid    = perimeter.centroid().cast<double>();
        const auto  to_centroid = [&](size_t k) { return (perimeter[k].cast<double>() - centroid).norm(); };
        double      nearest     = std::numeric_limits<double>::infinity();
        for (size_t k : candidates)
            nearest = std::min(nearest, to_centroid(k));
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [&](size_t k) { return to_centroid(k) > nearest + tolerance; }),
                         candidates.end());
    }
    size_t best = candidates.front();
    for (size_t k : candidates)
        if (rear_then_left(perimeter[k], perimeter[best]))
            best = k;
    return best;
}

// magnet apex: the vertex of the corner the modifier stands in, from its loop centroid `centre` and the
// nearest point `nearest`, or nullopt to keep that point. Q is the first crossing of a ray from the centre
// square to the line to P, within the limit as P is; the ray halfway between them picks which part of the
// loop from P to Q holds the corner; its vertex farthest from the centre wins when it lies deeper than P and Q.
// `modifier_distance` gives the distance of a point to the modifier's boundaries.
template<typename Distance>
static std::optional<MagnetPoint> find_apex_point(const Polygon &perimeter, const Lines &perimeter_lines, const Vec2d &centre,
                                                  const MagnetPoint &nearest, const Distance &modifier_distance,
                                                  std::optional<double> limit, double tolerance)
{
    const size_t count = perimeter.size();
    const Vec2d  to_p  = nearest.point.cast<double>() - centre;
    if (to_p.norm() <= tolerance)
        return std::nullopt;
    // The first crossing of a ray from the centre with the perimeter, beyond the tolerance.
    struct Hit {
        double            distance;
        Vec2d             point;
        PerimeterPosition position;
    };
    const auto cast = [&](const Vec2d &direction) {
        const Vec2d        u = direction.normalized();
        std::optional<Hit> best;
        for (size_t edge = 0; edge < count; ++edge) {
            const Vec2d  a     = perimeter_lines[edge].a.cast<double>();
            const Vec2d  ab    = perimeter_lines[edge].b.cast<double>() - a;
            const double denom = cross2(u, ab);
            if (denom == 0.) // Parallel to the ray.
                continue;
            const Vec2d  ao = a - centre;
            const double t  = cross2(ao, ab) / denom; // Along the ray.
            const double s  = cross2(ao, u) / denom;  // Along the edge.
            if (t <= tolerance || s < 0. || s > 1. || (best && t >= best->distance))
                continue;
            best = Hit{t, centre + t * u, s == 1. ? PerimeterPosition{(edge + 1) % count, 0.} : PerimeterPosition{edge, s}};
        }
        return best;
    };
    const Vec2d        side(-to_p.y(), to_p.x());
    std::optional<Hit> q = cast(side);
    if (const std::optional<Hit> other = cast(-side); other && (!q || other->distance < q->distance))
        q = other;
    if (!q || (limit && modifier_distance(q->point) > *limit + tolerance))
        return std::nullopt;
    const std::optional<Hit> r = cast(to_p.normalized() + (q->point - centre).normalized());
    if (!r)
        return std::nullopt;

    // The part from P to Q that holds R, taken forward: backward from P is forward from Q.
    const ArcLengths arcs(perimeter_lines);
    const double     sp    = arcs.along(nearest.position);
    const double     sq    = arcs.along(q->position);
    const bool       ahead = arcs.forward(sp, arcs.along(r->position)) <= arcs.forward(sp, sq);
    const double     ends  = std::max(to_p.norm(), q->distance);
    const std::optional<size_t> tip = ahead ? deepest_vertex(perimeter, arcs, sp, arcs.forward(sp, sq), centre, ends, tolerance) :
                                              deepest_vertex(perimeter, arcs, sq, arcs.forward(sq, sp), centre, ends, tolerance);
    if (!tip)
        return std::nullopt;
    return MagnetPoint{perimeter[*tip], PerimeterPosition{*tip, 0.}};
}

std::optional<MagnetPoint> find_touch_point(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                            const PerimeterSegment &touch, ModelVolumeType type, std::optional<double> limit,
                                            bool apex)
{
    if (!prepared.valid)
        return std::nullopt;
    const Polygon   &perimeter = prepared.polygon;
    const Lines      lines     = to_lines(perimeter);
    const ArcLengths arcs(lines);
    const double     tolerance = double(TOLERANCE_LINEAR);
    // The point the type takes on the crossing: its start for Left, its end for Right, its middle otherwise.
    // A segment runs forward from begin to end.
    const double begin  = arcs.along(touch.begin);
    double       middle = begin + arcs.forward(begin, arcs.along(touch.end)) / 2.;
    if (middle >= arcs.total())
        middle -= arcs.total();
    const Point       middle_point(arcs.point(lines, middle));
    PerimeterPosition position = type == ModelVolumeType::PRECISE_SEAM_LEFT  ? touch.begin :
                                 type == ModelVolumeType::PRECISE_SEAM_RIGHT ? touch.end :
                                                                               arcs.position(middle);
    if (position.parameter == 1.)
        position = {(position.edge_index + 1) % perimeter.size(), 0.};
    const Point seam = type == ModelVolumeType::PRECISE_SEAM_LEFT  ? touch.polyline.points.front() :
                       type == ModelVolumeType::PRECISE_SEAM_RIGHT ? touch.polyline.points.back() :
                                                                     middle_point;
    const MagnetPoint found{seam, position};
    if (!apex)
        return found;
    // apex looks from the centroid of the region that holds the touch, found by its middle: its ends lie on the
    // region's boundary.
    const auto region = std::find_if(modifier.begin(), modifier.end(), [&middle_point](const ModifierRegion &r) {
        return r.bounds.contains(middle_point) && r.polygon.contains(middle_point);
    });
    if (region == modifier.end() || region->polygon.contour.area() == 0.)
        return found;
    Lines modifier_lines;
    for (const ModifierRegion &r : modifier) {
        if (r.polygon.empty())
            continue;
        append(modifier_lines, to_lines(r.polygon.contour));
        for (const Polygon &hole : r.polygon.holes)
            append(modifier_lines, to_lines(hole));
    }
    const auto tree        = AABBTreeLines::build_aabb_tree_over_indexed_lines(modifier_lines);
    const auto to_modifier = [&](const Vec2d &point) {
        size_t       line = 0;
        Vec2d        nearest;
        const double squared = AABBTreeLines::squared_distance_to_indexed_lines(modifier_lines, tree, point, line, nearest);
        return squared < 0. ? std::numeric_limits<double>::infinity() : std::sqrt(squared);
    };
    return find_apex_point(perimeter, lines, region->polygon.contour.centroid().cast<double>(), found, to_modifier, limit, tolerance)
        .value_or(found);
}

std::optional<MagnetPoint> find_crossing_apex(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                              const SegmentExtraction &extracted)
{
    if (!prepared.valid || !extracted.valid || extracted.full_containment || extracted.discarded_fragments > 0 ||
        extracted.segments.size() < 2)
        return std::nullopt;
    const Polygon   &perimeter = prepared.polygon;
    const Lines      lines     = to_lines(perimeter);
    const ArcLengths arcs(lines);
    const double     tolerance = double(TOLERANCE_LINEAR);
    // The middle of each crossing, by the region it lies in. Segments run forward from begin to end.
    std::vector<std::vector<double>> middles(modifier.size());
    for (const PerimeterSegment &segment : extracted.segments) {
        const double begin  = arcs.along(segment.begin);
        double       middle = begin + arcs.forward(begin, arcs.along(segment.end)) / 2.;
        if (middle >= arcs.total())
            middle -= arcs.total();
        const Point point(arcs.point(lines, middle));
        const auto  region = std::find_if(modifier.begin(), modifier.end(), [&point](const ModifierRegion &r) {
            return r.bounds.contains(point) && r.polygon.contains(point);
        });
        // A crossing that no region holds cannot be placed: the type acts.
        if (region == modifier.end())
            return std::nullopt;
        middles[region - modifier.begin()].push_back(middle);
    }
    // A region crossed once makes the type act, apex or not.
    if (std::any_of(middles.begin(), middles.end(), [](const std::vector<double> &m) { return m.size() == 1; }))
        return std::nullopt;
    std::optional<MagnetPoint> best;
    double                     best_distance = std::numeric_limits<double>::infinity();
    for (size_t r = 0; r < modifier.size(); ++r) {
        std::vector<double> &m = middles[r];
        if (m.size() < 2 || modifier[r].polygon.contour.area() == 0.)
            continue;
        std::sort(m.begin(), m.end());
        // The shortest part that holds every middle leaves out the longest gap between neighbours.
        size_t gap = 0;
        for (size_t i = 1; i < m.size(); ++i)
            if (arcs.forward(m[i], m[(i + 1) % m.size()]) > arcs.forward(m[gap], m[(gap + 1) % m.size()]))
                gap = i;
        const double from   = m[(gap + 1) % m.size()];
        const double span   = arcs.total() - arcs.forward(m[gap], from);
        const Vec2d  centre = modifier[r].polygon.contour.centroid().cast<double>();
        const double ends   = std::max((arcs.point(lines, from) - centre).norm(), (arcs.point(lines, m[gap]) - centre).norm());
        const std::optional<size_t> tip = deepest_vertex(perimeter, arcs, from, span, centre, ends, tolerance);
        if (!tip)
            continue;
        // Tips of several regions: the one nearest the centroid of its own region, then the rear and left one.
        const double distance = (perimeter[*tip].cast<double>() - centre).norm();
        if (!best || distance < best_distance - tolerance ||
            (distance <= best_distance + tolerance && rear_then_left(perimeter[*tip], best->point))) {
            best          = MagnetPoint{perimeter[*tip], PerimeterPosition{*tip, 0.}};
            best_distance = distance;
        }
    }
    return best;
}

std::optional<MagnetPoint> find_magnet_point(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                             std::optional<double> limit, bool apex)
{
    if (!prepared.valid)
        return std::nullopt;
    const Polygon &perimeter = prepared.polygon;
    const size_t   count     = perimeter.size();
    static constexpr double infinity = std::numeric_limits<double>::infinity(); // Static: lambdas use it uncaptured.
    // Distances within the snapping radius count as equal. The limit takes the same margin, so that a touch
    // on an inclined edge, which computes slightly above zero, still counts for magnet0.
    const double tolerance = double(TOLERANCE_LINEAR);

    // Loops are the contours and holes of the regions in reach; each line remembers its loop.
    std::vector<const Polygon*> loops;
    for (const ModifierRegion &region : modifier) {
        if (region.polygon.empty())
            continue;
        if (limit) {
            const BoundingBox &a = prepared.bounds;
            const BoundingBox &b = region.bounds;
            // The boxes are at least as far apart as their gap along either axis.
            const double gap_x = std::max({0., double(b.min.x() - a.max.x()), double(a.min.x() - b.max.x())});
            const double gap_y = std::max({0., double(b.min.y() - a.max.y()), double(a.min.y() - b.max.y())});
            if (std::max(gap_x, gap_y) > *limit + tolerance)
                continue;
        }
        loops.push_back(&region.polygon.contour);
        for (const Polygon &hole : region.polygon.holes)
            loops.push_back(&hole);
    }
    Lines               modifier_lines;
    std::vector<size_t> line_loop;
    for (size_t loop = 0; loop < loops.size(); ++loop)
        for (const Line &line : to_lines(*loops[loop])) {
            modifier_lines.push_back(line);
            line_loop.push_back(loop);
        }
    if (modifier_lines.empty())
        return std::nullopt;
    const Lines perimeter_lines = to_lines(perimeter); // Edge i runs from vertex i to vertex i + 1.
    const auto  modifier_tree   = AABBTreeLines::build_aabb_tree_over_indexed_lines(modifier_lines);
    const auto  perimeter_tree  = AABBTreeLines::build_aabb_tree_over_indexed_lines(perimeter_lines);
    const auto  distance_to     = [](const Lines &lines, const auto &tree, const Vec2d &point) {
        size_t      line = 0;
        Vec2d       nearest;
        const double squared = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, point, line, nearest);
        return squared < 0. ? infinity : std::sqrt(squared);
    };

    // Stage 1: the nearest pair has a vertex of one contour, so vertex queries give the minimum distance.
    std::vector<double> vertex_distance(count);
    double              nearest = infinity;
    for (size_t k = 0; k < count; ++k) {
        vertex_distance[k] = distance_to(modifier_lines, modifier_tree, perimeter[k].cast<double>());
        nearest            = std::min(nearest, vertex_distance[k]);
    }
    struct ModifierVertex {
        Vec2d  point;
        size_t loop;
        double distance;
    };
    std::vector<ModifierVertex> modifier_vertices;
    for (size_t loop = 0; loop < loops.size(); ++loop)
        for (const Point &point : loops[loop]->points) {
            const Vec2d  p        = point.cast<double>();
            const double distance = distance_to(perimeter_lines, perimeter_tree, p);
            modifier_vertices.push_back({p, loop, distance});
            nearest = std::min(nearest, distance);
        }
    if (nearest == infinity || (limit && nearest > *limit + tolerance))
        return std::nullopt;
    const double threshold = std::min(nearest, limit ? *limit : infinity) + tolerance;
    // The tree truncates a query point to whole units on both axes, which moves it by up to sqrt(2): queries
    // take a radius 2 units larger, and every hit is then checked exactly against the threshold.
    const double search = (threshold + 2.) * (threshold + 2.);

    // Stage 2: every point of the minimum, each with the loops it is that near to.
    struct Candidate {
        Vec2d               point;
        size_t              edge;
        double              parameter;
        std::vector<size_t> loops;
    };
    std::vector<Candidate> candidates;
    const auto modifier_line_distance = [&](size_t line, const Vec2d &point) {
        return (project_on_segment(point, modifier_lines[line].a.cast<double>(), modifier_lines[line].b.cast<double>()).first - point).norm();
    };
    // (i) Perimeter vertices, with every loop within the threshold, not only the tree's single answer.
    for (size_t k = 0; k < count; ++k)
        if (vertex_distance[k] <= threshold) {
            const Vec2d         point = perimeter[k].cast<double>();
            std::vector<size_t> near;
            for (size_t line : AABBTreeLines::all_lines_in_radius(modifier_lines, modifier_tree, point, search))
                if (modifier_line_distance(line, point) <= threshold)
                    near.push_back(line_loop[line]);
            std::sort(near.begin(), near.end());
            near.erase(std::unique(near.begin(), near.end()), near.end());
            candidates.push_back({point, k, 0., std::move(near)});
        }
    // (ii) Projections of modifier vertices on every perimeter edge within the threshold.
    for (const ModifierVertex &vertex : modifier_vertices)
        if (vertex.distance <= threshold)
            for (size_t edge : AABBTreeLines::all_lines_in_radius(perimeter_lines, perimeter_tree, vertex.point, search)) {
                const auto [point, parameter] = project_on_segment(vertex.point, perimeter_lines[edge].a.cast<double>(),
                                                                   perimeter_lines[edge].b.cast<double>());
                if ((point - vertex.point).norm() <= threshold)
                    candidates.push_back({point, edge, parameter, {vertex.loop}});
            }
    // (iii) On parallel edges the minimum is an interval: add the projection of each loop centroid on the
    // edges of its candidates, where it is that near to the same loop.
    std::vector<std::optional<Vec2d>> centroids(loops.size());
    for (size_t loop = 0; loop < loops.size(); ++loop)
        if (loops[loop]->area() != 0.)
            centroids[loop] = loops[loop]->centroid().cast<double>();
    std::vector<std::pair<size_t, size_t>> edge_loops;
    for (const Candidate &candidate : candidates)
        for (size_t loop : candidate.loops) {
            edge_loops.emplace_back(candidate.edge, loop);
            // A candidate at a vertex also lies on the edge before it.
            if (candidate.parameter == 0.)
                edge_loops.emplace_back((candidate.edge + count - 1) % count, loop);
            else if (candidate.parameter == 1.)
                edge_loops.emplace_back((candidate.edge + 1) % count, loop);
        }
    std::sort(edge_loops.begin(), edge_loops.end());
    edge_loops.erase(std::unique(edge_loops.begin(), edge_loops.end()), edge_loops.end());
    for (const auto &[edge, loop] : edge_loops) {
        if (!centroids[loop])
            continue;
        const auto [point, parameter] = project_on_segment(*centroids[loop], perimeter_lines[edge].a.cast<double>(),
                                                           perimeter_lines[edge].b.cast<double>());
        // Measured to the lines of this loop only: another loop must not confirm a foreign centroid.
        const auto lines = AABBTreeLines::all_lines_in_radius(modifier_lines, modifier_tree, point, search);
        if (std::any_of(lines.begin(), lines.end(), [&, loop = loop, point = point](size_t line) {
                return line_loop[line] == loop && modifier_line_distance(line, point) <= threshold;
            }))
            candidates.push_back({point, edge, parameter, {loop}});
    }
    if (candidates.empty())
        return std::nullopt;

    // Each rule keeps the candidates within the tolerance of its best value, fixed before filtering.
    const auto keep_nearest = [&candidates, tolerance](const auto &value) {
        double best = infinity;
        for (const Candidate &candidate : candidates)
            best = std::min(best, value(candidate));
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                        [&](const Candidate &candidate) { return value(candidate) > best + tolerance; }),
                         candidates.end());
    };
    // The centroid of the candidate's own loop; skipped when some candidate has no loop with an area.
    const auto loop_centroid_distance = [&centroids](const Candidate &candidate) {
        double distance = infinity;
        for (size_t loop : candidate.loops)
            if (centroids[loop])
                distance = std::min(distance, (candidate.point - *centroids[loop]).norm());
        return distance;
    };
    if (std::all_of(candidates.begin(), candidates.end(),
                    [&](const Candidate &candidate) { return loop_centroid_distance(candidate) < infinity; }))
        keep_nearest(loop_centroid_distance);
    if (perimeter.area() != 0.) {
        const Vec2d centre = perimeter.centroid().cast<double>();
        keep_nearest([&centre](const Candidate &candidate) { return (candidate.point - centre).norm(); });
    }
    // Two different points always differ here; a repeat of one point keeps the first.
    const Candidate *best = nullptr;
    Point            best_point;
    for (const Candidate &candidate : candidates) {
        const Point point(candidate.point);
        if (best == nullptr || rear_then_left(point, best_point)) {
            best       = &candidate;
            best_point = point;
        }
    }
    const PerimeterPosition position = best->parameter == 1. ? PerimeterPosition{(best->edge + 1) % count, 0.} :
                                                               PerimeterPosition{best->edge, best->parameter};
    const MagnetPoint found{best_point, position};
    if (!apex)
        return found;
    // apex looks from the centroid of the modifier loop that the nearest point is nearest to.
    for (size_t loop : best->loops)
        if (centroids[loop]) {
            const auto to_modifier = [&](const Vec2d &point) { return distance_to(modifier_lines, modifier_tree, point); };
            return find_apex_point(perimeter, perimeter_lines, *centroids[loop], found, to_modifier, limit, tolerance).value_or(found);
        }
    return found;
}

// --- Drivers ---------------------------------------------------------------------------------------------------

// The type to apply on this layer, or nullopt where the rules exclude the layer. A modifier missing from
// `effective` keeps its own type on every layer.
static std::optional<ModelVolumeType> active_type(const EffectiveModifiers &effective, const ModelVolume *modifier, size_t layer_id)
{
    if (const auto it = effective.find(modifier); it != effective.end())
        return it->second.active_on(layer_id) ? std::optional(it->second.type) : std::nullopt;
    return modifier->type();
}

// The magnet of a modifier, or null.
static const EffectiveModifier *magnet_of(const EffectiveModifiers &effective, const ModelVolume *modifier)
{
    if (const auto it = effective.find(modifier); it != effective.end() && it->second.magnet)
        return &it->second;
    return nullptr;
}

// Whether a magnet acts instead of the type: the modifier does not cross the perimeter. Full containment
// comes first, as it keeps its closed segment; only discarded fragments still count as a crossing.
static bool magnet_acts(const SegmentExtraction &extracted)
{
    return extracted.valid && (extracted.full_containment || (extracted.segments.empty() && extracted.discarded_fragments == 0));
}

// The same for a magnet with apex, for which one crossing shorter than its length is only a touch: as the
// modifier nears a corner, one side is touched before the other, and the type must not take that moment.
static bool magnet_acts(const SegmentExtraction &extracted, const EffectiveModifier &magnet)
{
    return magnet_acts(extracted) ||
           (magnet.magnet_apex && extracted.valid && !extracted.full_containment && extracted.discarded_fragments == 0 &&
            extracted.segments.size() == 1 && extracted.segments.front().polyline.length() < magnet.magnet_apex_length);
}

// The point of a magnet where it acts: the nearest point, or for the touch of a magnet with apex the point its
// type takes on the touch.
static std::optional<MagnetPoint> acting_magnet_point(const PreparedPerimeter &prepared, const ModifierRegions &regions,
                                                      const SegmentExtraction &extracted, const EffectiveModifier &magnet,
                                                      ModelVolumeType type)
{
    if (magnet_acts(extracted))
        return find_magnet_point(prepared, regions, magnet.magnet_limit, magnet.magnet_apex);
    return find_touch_point(prepared, regions, extracted.segments.front(), type, magnet.magnet_limit, magnet.magnet_apex);
}

// Records one extraction for the warnings: the type of a modifier with a discarded fragment, the
// evaluated/reached flags for "had no effect" (a discarded fragment counts as reached, a contact not).
// `magnet` is set where a magnet acts instead of the type: then only a found point counts as reached.
static void record_extraction(PreciseSeamWarnings *warnings, const ModelVolume *modifier, ModelVolumeType type,
                              const SegmentExtraction &extracted, std::optional<bool> magnet = std::nullopt)
{
    if (warnings == nullptr || !extracted.valid)
        return;
    if (extracted.discarded_fragments > 0)
        PreciseSeamWarnings::mark(warnings->failed_types, type);
    const auto it = warnings->modifier_usage.find(modifier);
    if (it == warnings->modifier_usage.end())
        return;
    PreciseSeamWarnings::ModifierUsage &usage = it->second;
    PreciseSeamWarnings::set(usage.checked);
    if (magnet ? *magnet : (!extracted.segments.empty() || extracted.full_containment || extracted.discarded_fragments > 0))
        PreciseSeamWarnings::set(usage.reached);
}

std::optional<Point> insert_strong_seam_point_by_rules(const std::vector<const ModelVolume*> &strong_volumes, Polygon &polygon,
                                                       const PreparedPerimeter &prepared, const Layer *layer,
                                                       const ModifierRegionsCache &slices_cache, PreciseSeamWarnings *warnings,
                                                       const EffectiveModifiers &effective, bool quiet)
{
    // Mirrors insert_strong_seam_point() in PreciseSeam.cpp: keep the shared steps in step.
    assert(&prepared.polygon == &polygon);
    // A probe records nothing: the usage flags and warning masks belong to the real pass.
    assert(!quiet || warnings == nullptr);
    if (strong_volumes.empty() || layer == nullptr || !prepared.valid)
        return std::nullopt;
    const size_t layer_id = layer->id() - layer->object()->slicing_parameters().raft_layers();
    for (const ModelVolume *modifier : strong_volumes) {
        // Layer rules filter before any usage is recorded: an excluded layer is not "no effect".
        const std::optional<ModelVolumeType> active = active_type(effective, modifier, layer_id);
        if (!active)
            continue;
        const auto it = slices_cache.find(modifier);
        if (it == slices_cache.end() || layer_id >= it->second.size())
            continue;
        const ModelVolumeType type = *active;
        const SegmentExtraction extracted = extract_perimeter_segments(
            prepared, it->second[layer_id], type, {layer, modifier, warnings, quiet});
        // A magnet that does not cross the perimeter puts the seam at the nearest point, whatever the type.
        if (const EffectiveModifier *magnet = magnet_of(effective, modifier); magnet != nullptr && magnet_acts(extracted, *magnet)) {
            const std::optional<MagnetPoint> point = acting_magnet_point(prepared, it->second[layer_id], extracted, *magnet, type);
            record_extraction(warnings, modifier, type, extracted, point.has_value());
            if (!point)
                continue;
            const auto [seam_point, seam_index] = insert_point_into_perimeter(point->point, point->position.edge_index, polygon);
            refine_at_vertex(seam_index, +1, polygon);
            refine_at_vertex(seam_index, -1, polygon);
            return seam_point;
        }
        record_extraction(warnings, modifier, type, extracted);
        // apex: two or more crossings put the seam at the tip of the corner between them.
        if (const EffectiveModifier *magnet = magnet_of(effective, modifier); magnet != nullptr && magnet->magnet_apex)
            if (const std::optional<MagnetPoint> tip = find_crossing_apex(prepared, it->second[layer_id], extracted)) {
                const auto [seam_point, seam_index] = insert_point_into_perimeter(tip->point, tip->position.edge_index, polygon);
                refine_at_vertex(seam_index, +1, polygon);
                refine_at_vertex(seam_index, -1, polygon);
                return seam_point;
            }
        // Full containment leaves no intersection to place the point on: skipped with a warning.
        if (extracted.full_containment) {
            if (warnings)
                PreciseSeamWarnings::mark(warnings->full_containment, type);
            continue;
        }
        // No usable segment, even if all fragments were discarded: pass the turn to the next modifier.
        if (extracted.segments.empty())
            continue;
        if (warnings && extracted.segments.size() > 1)
            PreciseSeamWarnings::mark(warnings->multiple_intersections, type);

        // Lengths are compared exactly: for a symmetric modifier the chosen face may alternate between
        // layers, which is accepted (the user gets the multiple-intersections warning).
        StrongSeamTarget target{Point(0, 0), 0};
        double longest = -1.;
        for (const PerimeterSegment &segment : extracted.segments) {
            assert(segment.strong_target.has_value());
            const StrongSeamTarget &candidate = *segment.strong_target;
            const bool farther_or_left = rear_then_left(candidate.point, target.point);
            if (segment.length > longest || (segment.length == longest && farther_or_left)) {
                longest = segment.length;
                target = candidate;
            }
        }
        const auto [seam_point, seam_index] = insert_point_into_perimeter(target.point, target.edge_index, polygon);
        // Preserve the insertion order: refining before first would shift the seam index.
        refine_at_vertex(seam_index, +1, polygon);
        refine_at_vertex(seam_index, -1, polygon);
        return seam_point; // No later strong modifier or weak processing for this perimeter.
    }
    return std::nullopt;
}

std::vector<WeakModifierSegment> collect_weak_modifier_segments_by_rules(const std::vector<const ModelVolume*> &weak_volumes,
                                                                         Polygon &polygon, const PreparedPerimeter &prepared,
                                                                         const Layer *layer, const ModifierRegionsCache &slices_cache,
                                                                         PreciseSeamWarnings *warnings, const EffectiveModifiers &effective)
{
    // Mirrors collect_weak_modifier_segments() in PreciseSeam.cpp: keep the shared steps in step.
    assert(&prepared.polygon == &polygon);
    std::vector<WeakModifierSegment> result;
    if (weak_volumes.empty() || layer == nullptr || !prepared.valid)
        return result;

    // Modifier cache indices exclude raft layers, unlike Layer::id().
    const size_t raft_layers = layer->object()->slicing_parameters().raft_layers();
    const size_t layer_id = layer->id() - raft_layers;
    for (const ModelVolume *modifier_volume : weak_volumes) {
        const std::optional<ModelVolumeType> active = active_type(effective, modifier_volume, layer_id);
        if (!active)
            continue;
        const auto it = slices_cache.find(modifier_volume);
        if (it == slices_cache.end() || layer_id >= it->second.size())
            continue;
        const ModelVolumeType volume_type = *active;
        const ExtractionContext context{layer, modifier_volume, warnings};
        const SegmentExtraction extracted = extract_perimeter_segments(prepared, it->second[layer_id], volume_type, context);
        const auto type = convert_weak_modifier_type(volume_type);
        // A magnet that does not cross the perimeter types its nearest point only: a zone of one point.
        if (const EffectiveModifier *magnet = magnet_of(effective, modifier_volume); magnet != nullptr && magnet_acts(extracted, *magnet)) {
            const std::optional<MagnetPoint> point = acting_magnet_point(prepared, it->second[layer_id], extracted, *magnet, volume_type);
            record_extraction(warnings, modifier_volume, volume_type, extracted, point.has_value());
            if (point)
                result.push_back({type, point->point, point->position, point->point, point->position, false, modifier_volume});
            continue;
        }
        record_extraction(warnings, modifier_volume, volume_type, extracted);
        // apex: two or more crossings type the tip of the corner between them only.
        if (const EffectiveModifier *magnet = magnet_of(effective, modifier_volume); magnet != nullptr && magnet->magnet_apex)
            if (const std::optional<MagnetPoint> tip = find_crossing_apex(prepared, it->second[layer_id], extracted)) {
                result.push_back({type, tip->point, tip->position, tip->point, tip->position, false, modifier_volume});
                continue;
            }
        // Full containment works like seam painting all round: the whole perimeter takes the type.
        // A perimeter left whole Blocked also gets a warning (see apply_weak_modifiers_to_perimeter()).
        if (extracted.full_containment) {
            result.push_back({type, Point(), PerimeterPosition{0, 0.}, Point(), PerimeterPosition{0, 0.}, true, modifier_volume});
            continue;
        }
        for (const PerimeterSegment &segment : extracted.segments) {
            result.push_back({type, segment.polyline.points.front(), segment.begin,
                                   segment.polyline.points.back(), segment.end, false, modifier_volume});
        }
    }
    // Boundary insertion, refinement and last-write-wins priority use the established path.
    return prepare_weak_modifier_segments(std::move(result), polygon);
}

} // namespace PreciseSeam
} // namespace Slic3r
