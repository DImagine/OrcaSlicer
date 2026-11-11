#include "OuterWallPasses.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AABBTreeLines.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Flow.hpp"
#include "Line.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "Polyline.hpp"
#include "VariableWidth.hpp"
#include "libslic3r.h"

namespace Slic3r {

namespace {

using OutlineDistancer = AABBTreeLines::LinesDistancer<Line>;

// Outline to compare: either as sliced, or with the holes filled.
ExPolygons comparable_outline(const ExPolygons &src, bool with_holes)
{
    if (with_holes)
        return src;
    Polygons contours;
    contours.reserve(src.size());
    for (const ExPolygon &expoly : src)
        contours.emplace_back(expoly.contour);
    return union_ex(contours);
}

// Distance of a point from an outline, positive inside. Returns false if the outline is empty.
bool signed_distance(const OutlineDistancer &outline, const Point &point, double &distance, Vec2d &nearest)
{
    const auto [dist, line_idx, nearest_point] = outline.distance_from_lines_extra<false>(point);
    if (! std::isfinite(dist))
        return false;
    distance = outline.outside(point) > 0 ? - dist : dist;
    nearest  = nearest_point;
    return true;
}

// The ring, a contour or a hole, of every line of an outline, in the order of to_lines().
std::vector<size_t> rings_of_lines(const ExPolygons &outline)
{
    std::vector<size_t> rings;
    size_t              ring = 0;
    for (const ExPolygon &expoly : outline)
        for (size_t i = 0; i <= expoly.holes.size(); ++ i, ++ ring)
            rings.insert(rings.end(), (i == 0 ? expoly.contour : expoly.holes[i - 1]).points.size(), ring);
    return rings;
}

// Where the material of an outline ends on the way out of it from a point, in the units of the outward vector,
// between the given two. Finds the edge a loop follows where another edge of the same ring is nearer to the point:
// the material has to end there once only, on the given ring. Returns false otherwise.
bool edge_on_way_out(const OutlineDistancer &outline, const std::vector<size_t> &rings, size_t ring, const Vec2d &point,
                     const Vec2d &outward, double lowest, double highest, double &at)
{
    const Line way(Point((point + lowest * outward).cast<coord_t>()), Point((point + highest * outward).cast<coord_t>()));
    // Crossings of the way with the outline.
    struct Crossing
    {
        size_t ring;
        double at;
        // Whether the material ends there, and whether the way goes through an end of the edge.
        bool   out;
        bool   vertex;
        bool operator<(const Crossing &rhs) const { return this->ring < rhs.ring || (this->ring == rhs.ring && this->at < rhs.at); }
    };
    const double          same_point = scaled<double>(0.0001);
    std::vector<Crossing> crossings;
    for (const auto &[crossing, line_idx] : outline.intersections_with_line<false>(way)) {
        // The material lies to the left of an edge of the outlines: it ends where the way out leaves to the right.
        const Line  &edge = outline.get_line(line_idx);
        const Vec2d  dir  = (edge.b - edge.a).cast<double>();
        const double turn = dir.x() * outward.y() - dir.y() * outward.x();
        if (turn == 0.)
            continue;
        const Vec2d where = crossing.cast<double>();
        crossings.push_back({ rings[line_idx], (where - point).dot(outward) / outward.squaredNorm(), turn < 0.,
                              (where - edge.a.cast<double>()).norm() <= same_point || (where - edge.b.cast<double>()).norm() <= same_point });
    }
    std::sort(crossings.begin(), crossings.end());
    size_t ways_out = 0;
    for (size_t begin = 0; begin < crossings.size();) {
        // Through a vertex the way crosses both its edges at one point. Where it only touches the vertex, it leaves
        // the material and enters it at once: that is no way out.
        size_t end  = begin;
        int    ends = 0;
        for (; end < crossings.size() && crossings[end].ring == crossings[begin].ring &&
               (crossings[end].at - crossings[begin].at) * outward.norm() <= same_point; ++ end)
            ends += crossings[end].out ? 1 : -1;
        // A vertex crossed along one of its edges only: the other one runs along the way, which tells nothing.
        if (end == begin + 1 && crossings[begin].vertex)
            return false;
        if (ends > 0) {
            // A hole, a notch or another part beside the wall gives the way more than one end of the material, or
            // one on another ring: which of them the loop follows is not known.
            if (++ ways_out > 1 || crossings[begin].ring != ring)
                return false;
            at = crossings[begin].at;
        }
        begin = end;
    }
    return ways_out == 1;
}

// Point of a wall loop the passes are derived from.
struct LoopSample
{
    Point  point;
    // Unit vector pointing out of the material, scaled so that moving the point along it by d offsets the adjacent
    // edges by d (miter offset of a polygon vertex).
    Vec2d  outward;
    // Offset of each pass outline from the middle outline at this point, positive outwards, scaled.
    std::vector<double> offsets;
    // Shift of the middle line of each pass of this wall outwards and what its width gains, scaled.
    std::vector<double> shifts;
    std::vector<double> gains;
    // Whether this is a vertex of the loop or a point added inside an edge.
    bool   vertex;
    // Whether the wall follows the surface of the object here.
    bool   on_surface { false };
    // Whether the wall generator took the wall for hanging in the air here.
    bool   overhang { false };
};

// How much narrower and wider than its regular width a pass of a wall may get, scaled: within the range the walls keep
// to as long as they can, and within the limits.
struct WallReach
{
    double inwards { 0. };
    double outwards { 0. };
    double limit_inwards { 0. };
    double limit_outwards { 0. };
};

// What the walls of a layer region printed in passes share, the outer wall first.
struct WallSet
{
    // Regular width of a wall, unscaled.
    std::vector<double>    widths;
    // Distance of the middle line of a wall from the surface of the object, scaled.
    std::vector<double>    depths;
    std::vector<WallReach> reach;
    // Flow of a wall at the height of a pass.
    std::vector<Flow>      pass_flows;
    // Where the middle line of an inner wall runs along the surface of the object: the outline of the layer offset
    // by the depth of the wall. Empty for the outer wall.
    std::vector<OutlineDistancer> surface_lines;

    size_t count() const { return this->widths.size(); }
};

// Loop of a wall taken for printing in passes.
struct WallLoop
{
    // Where the loop is stored among the perimeters.
    ExtrusionEntity           **slot { nullptr };
    // Number of the wall, zero for the outer one.
    size_t                      wall { 0 };
    Polygon                     polygon;
    // 1 if the material is to the left of the loop, -1 if to the right.
    double                      side { 1. };
    // Whether the loop runs around a hole of the material.
    bool                        hole { false };
    std::vector<LoopSample>     samples;
    std::vector<ExtrusionLoop>  passes;
    bool                        valid { false };
    // TEMP diagnostics: why the loop is left as it is.
    std::string                 reject;

    const ExtrusionLoop& loop() const { return *static_cast<const ExtrusionLoop*>(*this->slot); }
};

// Sample a wall loop and measure the offsets of the pass outlines along it. Returns false if the loop has to be
// printed the regular way.
bool sample_loop(WallLoop                            &wall_loop,
                 const WallSet                       &walls,
                 const OuterWallPassSlices           &slices,
                 const OutlineDistancer              &layer_outline,
                 const OutlineDistancer              &mid_outline,
                 const std::vector<OutlineDistancer> &pass_outlines,
                 const std::vector<size_t>           &mid_rings,
                 const std::vector<std::vector<size_t>> &pass_rings,
                 bool                                 with_holes)
{
    const ExtrusionLoop &loop          = wall_loop.loop();
    const size_t         wall          = wall_loop.wall;
    auto                 refuse        = [&wall_loop](const std::string &why) { wall_loop.reject = why; return false; };
    const double         regular_width = walls.widths[wall];
    const ExtrusionRole  role          = wall == 0 ? erExternalPerimeter : erPerimeter;
    // Only a wall of its regular width is split, and a wall held back until the infill is printed must not be
    // printed first. The stretches the wall generator took for hanging in the air are split with the rest of the
    // loop: the wall is taken for hanging by the whole layer, while a pass reaches out half as far over the pass
    // below it. They keep their role, that is their speed and cooling, but not their flow.
    if (loop.print_after_infill)
        return refuse("after_infill");
    Lines overhangs;
    for (const ExtrusionPath &path : loop.paths) {
        if (path.role() == erOverhangPerimeter) {
            append(overhangs, path.polyline.to_polyline().lines());
            continue;
        }
        if (path.role() != role || std::abs(double(path.width) - regular_width) > 1e-3)
            return refuse("width_or_role");
    }
    const OutlineDistancer overhang_lines(std::move(overhangs));

    wall_loop.polygon = loop.polygon();
    wall_loop.polygon.remove_duplicate_points();
    const Points &points = wall_loop.polygon.points;
    if (points.size() < 3)
        return refuse("degenerate");

    // Sample the loop densely enough to follow a slope changing along a straight edge.
    const double             max_segment = scaled<double>(0.5);
    const size_t             num_passes  = slices.count();
    std::vector<LoopSample> &samples     = wall_loop.samples;
    std::vector<Vec2d>       edge_normals;
    samples.clear();
    edge_normals.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++ i) {
        const Vec2d dir = (points[(i + 1) % points.size()] - points[i]).cast<double>().normalized();
        // Normal to the right of the edge.
        edge_normals.emplace_back(dir.y(), - dir.x());
    }
    for (size_t i = 0; i < points.size(); ++ i) {
        const Vec2d  a      = points[i].cast<double>();
        const Vec2d  ab     = points[(i + 1) % points.size()].cast<double>() - a;
        const size_t steps  = std::max<size_t>(1, size_t(std::ceil(ab.norm() / max_segment)));
        const Vec2d &n_next = edge_normals[i];
        const Vec2d &n_prev = edge_normals[(i + points.size() - 1) % points.size()];
        for (size_t j = 0; j < steps; ++ j) {
            LoopSample sample;
            sample.vertex = j == 0;
            sample.point  = sample.vertex ? points[i] : Point((a + ab * (double(j) / double(steps))).cast<coord_t>());
            if (sample.vertex) {
                // Miter offset, limited at sharp corners.
                const double denom = 1. + n_prev.dot(n_next);
                sample.outward = denom > 0.25 ? Vec2d((n_prev + n_next) / denom) : Vec2d(2. * (n_prev + n_next).normalized());
                if (! std::isfinite(sample.outward.x()) || ! std::isfinite(sample.outward.y()))
                    sample.outward = n_next;
            } else
                sample.outward = n_next;
            {
                // The segment from this point to the next one hangs if its middle lies on a hanging path: the ends
                // of a hanging path are the ends of the regular paths beside it as well.
                const Point middle((a + ab * ((double(j) + 0.5) / double(steps))).cast<coord_t>());
                const auto [dist, line_idx, nearest_point] = overhang_lines.distance_from_lines_extra<false>(middle);
                sample.overhang = std::isfinite(dist) && dist <= scaled<double>(0.01);
            }
            samples.emplace_back(std::move(sample));
        }
    }

    // Measure the offsets of the pass outlines and find out, which side of the loop the material is on.
    // The material is inside a contour and outside a hole. The normals point to the right of the loop.
    const bool is_hole = (loop.loop_role() & elrHole) == elrHole;
    if (! with_holes && is_hole)
        return refuse("holes_off");
    const double side = wall_loop.polygon.is_counter_clockwise() != is_hole ? 1. : -1.;
    wall_loop.side = side;
    wall_loop.hole = is_hole;
    // The offsets are measured against the nearest edge of the outlines, which has to be the edge the loop follows:
    // the distance from the middle outline, positive inside, has to fall towards the outer side of the loop.
    // It is probed on both sides of the loop, which works even for a loop running along the outline itself, as it
    // does with an XY size compensation of half the wall width. A loop the compensations moved beside another edge
    // of the outlines, such as the hole of a thin wall they thickened, is left as it is.
    const double side_probe     = scaled<double>(0.05);
    bool         side_confirmed = false;
    // How far all the walls together move the outer edge at most, for the points where it is looked up on the way
    // out of the material.
    const double edge_tolerance = scaled<double>(0.05);
    double       most_inwards   = edge_tolerance;
    double       most_outwards  = edge_tolerance;
    for (const WallReach &reach : walls.reach) {
        most_inwards  += reach.limit_inwards;
        most_outwards += reach.limit_outwards;
    }
    // Whether the distance from an outline falls towards the outer side of the loop at a point: 1 if it does,
    // -1 if it rises, 0 if the probes on both sides of the point do not tell.
    auto outline_side = [side, side_probe](const OutlineDistancer &outline, const Vec2d &point, const Vec2d &outward) {
        const Vec2d right = side_probe * outward.normalized();
        double      distance_right;
        double      distance_left;
        Vec2d       nearest;
        if (! signed_distance(outline, Point((point + right).cast<coord_t>()), distance_right, nearest) ||
            ! signed_distance(outline, Point((point - right).cast<coord_t>()), distance_left, nearest) ||
            std::abs(distance_left - distance_right) <= 0.5 * side_probe)
            return 0;
        return (distance_left > distance_right) == (side > 0.) ? 1 : -1;
    };
    // Where a wall follows the surface of the object, the outline of the layer lies at the depth of the wall off
    // its middle on its outer side, at a corner along the miter: the distance alone would take the corner of
    // a hole, further from the outline than the sides of the loop, for something else. Where the outline is not
    // there, the loop runs along the border of another region of the layer: it is neither offset nor checked
    // against the outlines, whose nearest edge is the far side of the object.
    const double surface_offset    = side * walls.depths[wall];
    const double surface_tolerance = scaled<double>(0.05);
    // An inner wall is measured where the outer wall runs in front of it, so that the walls sharing an offset
    // take the same one.
    const double measure_offset    = side * (walls.depths[wall] - walls.depths.front());
    for (LoopSample &sample : samples) {
        const Vec2d point = sample.point.cast<double>();
        double      surface_distance;
        Vec2d       surface_nearest;
        // The sides of an inner wall run past the ends of the sides of a hole it goes around, where a point moved
        // by the depth of the wall misses the outline. An inner wall is thus looked up on the outline offset by
        // its depth, as the wall generator makes it.
        if (wall == 0 ? ! signed_distance(layer_outline, Point((point + surface_offset * sample.outward).cast<coord_t>()), surface_distance,
                                          surface_nearest) :
                        ! signed_distance(walls.surface_lines[wall], sample.point, surface_distance, surface_nearest))
            return refuse("no_surface_distance");
        if (std::abs(surface_distance) > surface_tolerance) {
            sample.offsets.assign(num_passes, 0.);
            continue;
        }
        const Vec2d measured_at = wall == 0 ? point : Vec2d(point + measure_offset * sample.outward);
        const Point measured(measured_at.cast<coord_t>());
        double      mid_distance;
        Vec2d       mid_nearest;
        if (! signed_distance(mid_outline, wall == 0 ? sample.point : measured, mid_distance, mid_nearest))
            return refuse("no_mid_distance");
        const int mid_side = outline_side(mid_outline, measured_at, sample.outward);
        if (mid_side < 0)
            return refuse("mid_side");
        side_confirmed |= mid_side > 0;
        sample.on_surface = true;
        sample.offsets.reserve(num_passes);
        for (const OutlineDistancer &pass_outline : pass_outlines) {
            double pass_distance;
            Vec2d  pass_nearest;
            if (! signed_distance(pass_outline, wall == 0 ? sample.point : measured, pass_distance, pass_nearest))
                return refuse("no_pass_distance");
            // The nearest edge may be another one in the outline of a pass only. Where it is the far side of a part
            // tapering off, the same ring of the outline coming round, the edge the loop follows is looked up on
            // the way out of the material, in the middle outline as well, as far as the walls may move it. The
            // planning keeps the edge within that, so the only end of the material there is the one the loop
            // follows. A loop beside an edge of another ring, a hole or another part, or with none or several ends
            // of the material on the way out is left as it is: which edge it follows there is not known.
            if (outline_side(pass_outline, measured_at, sample.outward) < 0) {
                const Vec2d                way_out = side * sample.outward;
                const Point                from(measured_at.cast<coord_t>());
                const std::vector<size_t> &rings   = pass_rings[&pass_outline - pass_outlines.data()];
                double                     mid_at  = 0.;
                double                     pass_at = 0.;
                if (! edge_on_way_out(mid_outline, mid_rings, mid_rings[std::get<1>(mid_outline.distance_from_lines_extra<false>(from))],
                                      measured_at, way_out, 0., 4. * mid_distance + edge_tolerance, mid_at) ||
                    ! edge_on_way_out(pass_outline, rings, rings[std::get<1>(pass_outline.distance_from_lines_extra<false>(from))],
                                      measured_at, way_out, mid_at - most_inwards, mid_at + most_outwards, pass_at))
                    return refuse("pass_side");
                sample.offsets.emplace_back(pass_at - mid_at);
                continue;
            }
            sample.offsets.emplace_back(pass_distance - mid_distance);
        }
    }
    return side_confirmed || refuse("side_unconfirmed");
}

// Split the offsets of the pass outlines between the walls running one behind another at each point of the loop.
// wall_lines: for this wall, the loops of every other wall printed in passes moved to where this wall runs along them.
// Returns false if the loop takes none: an inner wall with no wall printed in passes in front of it, or one that
// has nothing to take, which is then better printed in one piece.
bool share_offsets(WallLoop &wall_loop, const WallSet &walls, const std::vector<OutlineDistancer> &wall_lines, size_t num_passes)
{
    const size_t wall      = wall_loop.wall;
    const double tolerance = scaled<double>(0.05);
    const double negligible = scaled<double>(0.001);
    bool         shares    = false;
    for (LoopSample &sample : wall_loop.samples) {
        sample.shifts.assign(num_passes, 0.);
        sample.gains.assign(num_passes, 0.);
        if (! sample.on_surface)
            continue;
        // Whether another wall printed in passes runs along this one here.
        auto wall_at = [&wall_lines, &sample, tolerance](size_t other) {
            const auto [dist, line_idx, nearest_point] = wall_lines[other].distance_from_lines_extra<false>(sample.point);
            return std::isfinite(dist) && dist <= tolerance;
        };
        // The walls in front of this one and those behind it, as far as they run here one behind another.
        bool fronted = true;
        for (size_t other = 0; other < wall && fronted; ++ other)
            fronted = wall_at(other);
        if (! fronted)
            continue;
        size_t last = wall;
        while (last + 1 < walls.count() && wall_at(last + 1))
            ++ last;
        for (size_t pass = 0; pass < num_passes; ++ pass) {
            const double offset  = sample.offsets[pass];
            const bool   out     = offset >= 0.;
            double       within  = 0.;
            double       limit   = 0.;
            for (size_t other = 0; other <= last; ++ other) {
                within += out ? walls.reach[other].outwards : walls.reach[other].inwards;
                limit  += out ? walls.reach[other].limit_outwards : walls.reach[other].limit_inwards;
            }
            // The walls reach the end of the range they keep to together, and then their limits together.
            const double amount = std::min(std::abs(offset), limit);
            auto share = [&](size_t other) {
                const double reach_within = out ? walls.reach[other].outwards : walls.reach[other].inwards;
                const double reach_limit  = out ? walls.reach[other].limit_outwards : walls.reach[other].limit_inwards;
                if (amount <= within)
                    return within > 0. ? amount * (reach_within / within) : 0.;
                return reach_within + (amount - within) * ((reach_limit - reach_within) / (limit - within));
            };
            // The wall moves by what the walls behind it take and by a half of what it takes itself.
            double behind = 0.;
            for (size_t other = wall + 1; other <= last; ++ other)
                behind += share(other);
            const double own = share(wall);
            sample.shifts[pass] = (out ? 1. : -1.) * (behind + 0.5 * own);
            sample.gains[pass]  = (out ? 1. : -1.) * own;
            shares |= own > negligible;
        }
    }
    return shares || wall == 0;
}

// Build the passes of a wall loop from the shares of its samples. Returns false if a pass makes no loop.
bool build_passes(WallLoop &wall_loop, const WallSet &walls, size_t num_passes)
{
    const ExtrusionLoop           &loop          = wall_loop.loop();
    const size_t                   wall          = wall_loop.wall;
    const Flow                    &pass_flow     = walls.pass_flows[wall];
    const std::vector<LoopSample> &samples       = wall_loop.samples;
    const double                   regular_width = scaled<double>(walls.widths[wall]);
    // thick_polyline_to_multi_path() takes the width of the rectangle of the same cross section as the extrusion.
    const double width_correction = scaled<double>(pass_flow.height() * (1. - 0.25 * PI));
    const double linear_tolerance = scaled<double>(0.005);
    const float  width_tolerance  = scaled<float>(0.02);
    wall_loop.passes.clear();
    wall_loop.passes.reserve(num_passes);
    const ExtrusionRole role = wall == 0 ? erExternalPerimeter : erPerimeter;
    for (size_t pass = 0; pass < num_passes; ++ pass) {
        ThickPolyline polyline;
        std::vector<coordf_t> widths;
        std::vector<char>     hanging;
        auto append = [&](const LoopSample &sample) {
            polyline.points.emplace_back((sample.point.cast<double>() + wall_loop.side * sample.shifts[pass] * sample.outward).cast<coord_t>());
            widths.emplace_back(regular_width + sample.gains[pass] - width_correction);
            hanging.emplace_back(sample.overhang);
        };
        for (size_t i = 0; i < samples.size();) {
            // Samples of one edge: [i, edge_end), the next vertex is at edge_end.
            size_t edge_end = i + 1;
            while (edge_end < samples.size() && ! samples[edge_end].vertex)
                ++ edge_end;
            append(samples[i]);
            // Keep the points inside the edge only if the pass does not change linearly along it.
            const LoopSample &sample_a = samples[i];
            const LoopSample &sample_b = samples[edge_end % samples.size()];
            bool              linear   = true;
            for (size_t j = i + 1; j < edge_end && linear; ++ j) {
                const double t = double(j - i) / double(edge_end - i);
                // An overhanging stretch may begin and end inside the edge.
                linear = samples[j].overhang == sample_a.overhang &&
                         std::abs(samples[j].shifts[pass] - (sample_a.shifts[pass] + t * (sample_b.shifts[pass] - sample_a.shifts[pass]))) <= 0.5 * linear_tolerance &&
                         std::abs(samples[j].gains[pass] - (sample_a.gains[pass] + t * (sample_b.gains[pass] - sample_a.gains[pass]))) <= linear_tolerance;
            }
            if (! linear)
                for (size_t j = i + 1; j < edge_end; ++ j)
                    append(samples[j]);
            i = edge_end;
        }
        // Close the loop.
        polyline.points.emplace_back(polyline.points.front());
        widths.emplace_back(widths.front());
        hanging.emplace_back(hanging.front());

        // Split the pass into paths of a constant width, the stretches taken for hanging in the air apart from the
        // rest: a point tells whether the segment beginning at it hangs. The steps of the width are finer than
        // those of thin walls to keep the outer surface smooth.
        ExtrusionPaths paths;
        for (size_t begin = 0; begin + 1 < polyline.points.size();) {
            const bool hangs = hanging[begin];
            size_t     end   = begin + 1;
            while (end + 1 < polyline.points.size() && hanging[end] == hangs)
                ++ end;
            ThickPolyline stretch;
            stretch.points.assign(polyline.points.begin() + begin, polyline.points.begin() + end + 1);
            stretch.width.reserve(2 * (end - begin));
            for (size_t i = begin; i < end; ++ i) {
                stretch.width.emplace_back(widths[i]);
                stretch.width.emplace_back(widths[i + 1]);
            }
            Slic3r::append(paths, thick_polyline_to_multi_path(stretch, hangs ? erOverhangPerimeter : role, pass_flow, width_tolerance,
                                                       width_tolerance).paths);
            begin = end;
        }
        if (paths.empty() || paths.front().first_point() != paths.back().last_point())
            return false;
        ExtrusionLoop pass_loop(std::move(paths), loop.loop_role());
        pass_loop.inset_idx       = loop.inset_idx;
        pass_loop.outer_wall_pass = true;
        wall_loop.passes.emplace_back(std::move(pass_loop));
    }
    return true;
}

// The loops of the walls which may be printed in passes, in the order they are stored in.
void collect_wall_loops(ExtrusionEntityCollection &collection, size_t num_walls, std::vector<WallLoop> &wall_loops)
{
    for (ExtrusionEntity *&entity : collection.entities) {
        if (entity->is_collection()) {
            collect_wall_loops(*static_cast<ExtrusionEntityCollection*>(entity), num_walls, wall_loops);
        } else if (entity->is_loop() && entity->inset_idx >= 0 && size_t(entity->inset_idx) < num_walls) {
            WallLoop wall_loop;
            wall_loop.slot = &entity;
            wall_loop.wall = size_t(entity->inset_idx);
            wall_loops.emplace_back(std::move(wall_loop));
        }
    }
}

OuterWallPassReach reach_of_range(double min_width, double max_width)
{
    const double nominal_width = 0.5 * (min_width + max_width);
    return { nominal_width - min_width, max_width - nominal_width };
}

} // namespace

// Slices taken at different heights never match exactly: the contours are simplified with this precision.
static const float outline_tolerance = scaled<float>(0.003);
// Slivers the polygon clipping leaves along edges that coincide. Kept far below the tolerance, as whatever
// the opening removes adds to it.
static const float outline_sliver    = scaled<float>(0.0005);

bool outer_wall_passes_fit(const ExPolygons              &mid,
                           const std::vector<ExPolygons> &passes,
                           double                         min_width,
                           double                         max_width,
                           bool                           with_holes)
{
    return outer_wall_passes_fit(mid, passes, reach_of_range(min_width, max_width), with_holes);
}

bool outer_wall_passes_fit(const ExPolygons              &mid,
                           const std::vector<ExPolygons> &passes,
                           const OuterWallPassReach      &reach,
                           bool                           with_holes)
{
    const ExPolygons mid_outline = comparable_outline(mid, with_holes);
    // The widest pass may reach this far out and the narrowest pass has to reach at least this far.
    const Polygons   outer_limit = offset(mid_outline, scaled<float>(reach.outwards) + outline_tolerance);
    const ExPolygons inner_limit = offset_ex(mid_outline, - scaled<float>(reach.inwards) - outline_tolerance);
    for (const ExPolygons &pass : passes) {
        const ExPolygons pass_outline = comparable_outline(pass, with_holes);
        if (! opening_ex(diff_ex(pass_outline, outer_limit), outline_sliver).empty() ||
            ! opening_ex(diff_ex(inner_limit, pass_outline), outline_sliver).empty())
            return false;
    }
    return true;
}

bool outer_wall_pass_steps_fit(const std::vector<ExPolygons> &passes, const ExPolygons *above, double max_step, bool with_holes)
{
    // What is left of an outline shrunk by the step has to lie under the outline above it.
    auto covered = [max_step, with_holes](const ExPolygons &lower, const ExPolygons &upper) {
        const ExPolygons kept = offset_ex(comparable_outline(lower, with_holes), - scaled<float>(max_step) - outline_tolerance);
        return opening_ex(diff_ex(kept, comparable_outline(upper, with_holes)), outline_sliver).empty();
    };
    for (size_t pass = 0; pass + 1 < passes.size(); ++ pass)
        if (! covered(passes[pass], passes[pass + 1]))
            return false;
    return above == nullptr || passes.empty() || covered(passes.back(), *above);
}

void split_outer_wall_loops(ExtrusionEntityCollection              &perimeters,
                            std::vector<ExtrusionEntityCollection> &lower_passes,
                            const OuterWallPassSlices              &slices,
                            const ExPolygons                       &layer_outline,
                            const Flow                             &pass_flow,
                            bool                                    with_holes,
                            const OuterWallPassInnerWalls          *inner_walls,
                            std::vector<std::string>               *rejections)
{
    const size_t num_passes = slices.count();
    if (num_passes < 2)
        return;
    lower_passes.resize(num_passes - 1);
    const OutlineDistancer        layer_outline_lines(to_lines(layer_outline));
    const OutlineDistancer        mid_outline(to_lines(slices.mid));
    std::vector<OutlineDistancer> pass_outlines;
    pass_outlines.reserve(num_passes);
    for (const ExPolygons &pass : slices.passes)
        pass_outlines.emplace_back(to_lines(pass));
    // The rings of the lines of these outlines.
    const std::vector<size_t>        mid_rings = rings_of_lines(slices.mid);
    std::vector<std::vector<size_t>> pass_rings;
    pass_rings.reserve(num_passes);
    for (const ExPolygons &pass : slices.passes)
        pass_rings.emplace_back(rings_of_lines(pass));

    // The walls printed in passes: the outer wall of the nominal width and the inner walls behind it.
    WallSet      walls;
    const size_t num_walls = inner_walls != nullptr ? std::max<size_t>(1, slices.walls) : 1;
    for (size_t wall = 0; wall < num_walls; ++ wall) {
        const double width = wall == 0 ? slices.nominal_width() : double(inner_walls->pass_flow.width());
        walls.widths.emplace_back(width);
        walls.depths.emplace_back(wall == 0 ? 0.5 * scaled<double>(width) :
                                  walls.depths.back() + scaled<double>(wall == 1 ? inner_walls->first_spacing : inner_walls->spacing));
        walls.pass_flows.emplace_back(wall == 0 ? pass_flow : inner_walls->pass_flow);
        WallReach reach;
        reach.inwards        = scaled<double>(std::max(0., width - slices.min_width));
        reach.outwards       = scaled<double>(std::max(0., slices.max_width - width));
        reach.limit_inwards  = std::max(reach.inwards, scaled<double>(width - slices.lowest_width()));
        reach.limit_outwards = std::max(reach.outwards, scaled<double>(slices.highest_width() - width));
        walls.reach.emplace_back(reach);
        walls.surface_lines.emplace_back(wall == 0 ? Lines() : to_lines(offset_ex(layer_outline, - float(walls.depths.back()))));
    }

    std::vector<WallLoop> wall_loops;
    collect_wall_loops(perimeters, num_walls, wall_loops);
    for (WallLoop &wall_loop : wall_loops)
        wall_loop.valid = sample_loop(wall_loop, walls, slices, layer_outline_lines, mid_outline, pass_outlines, mid_rings, pass_rings,
                                      with_holes);
    // The walls share the offsets only with the walls printed in passes: once a loop drops out, the shares of the
    // loops beside it change.
    for (bool changed = true; changed;) {
        changed = false;
        // For every wall, the loops of every other wall printed in passes moved to where the wall runs along them:
        // offset by the difference of their depths, as the wall generator offsets one wall from another. Moving
        // a point of a wall along its normal would miss the other wall beside a corner of a hole, where the sides
        // of the walls differ in length.
        std::vector<std::vector<OutlineDistancer>> wall_lines(num_walls);
        for (size_t wall = 0; wall < num_walls; ++ wall) {
            wall_lines[wall].reserve(num_walls);
            for (size_t other = 0; other < num_walls; ++ other) {
                Lines lines;
                if (other != wall)
                    for (const WallLoop &wall_loop : wall_loops)
                        if (wall_loop.valid && wall_loop.wall == other) {
                            // Deeper into the material a contour shrinks and a hole grows.
                            Polygon polygon = wall_loop.polygon;
                            polygon.make_counter_clockwise();
                            const double deeper = walls.depths[wall] - walls.depths[other];
                            append(lines, to_lines(offset(polygon, float(wall_loop.hole ? deeper : - deeper))));
                        }
                wall_lines[wall].emplace_back(std::move(lines));
            }
        }
        for (WallLoop &wall_loop : wall_loops)
            if (wall_loop.valid && ! share_offsets(wall_loop, walls, wall_lines[wall_loop.wall], num_passes)) {
                wall_loop.reject = "no_share";
                wall_loop.valid = false;
                changed         = true;
            }
        if (changed)
            continue;
        for (WallLoop &wall_loop : wall_loops)
            if (wall_loop.valid && ! build_passes(wall_loop, walls, num_passes)) {
                wall_loop.reject = "build";
                wall_loop.valid = false;
                changed         = true;
            }
    }

    for (WallLoop &wall_loop : wall_loops) {
        if (! wall_loop.valid) {
            if (rejections != nullptr && wall_loop.wall == 0) {
                const Point at = wall_loop.loop().first_point();
                rejections->emplace_back((wall_loop.reject.empty() ? std::string("unknown") : wall_loop.reject) +
                    " len " + std::to_string(unscaled<double>(wall_loop.loop().length())) +
                    " at " + std::to_string(unscaled<double>(at.x())) + " " + std::to_string(unscaled<double>(at.y())));
            }
            continue;
        }
        for (size_t pass = 0; pass + 1 < num_passes; ++ pass)
            lower_passes[pass].entities.emplace_back(new ExtrusionLoop(std::move(wall_loop.passes[pass])));
        // The top pass is printed at the layer height and takes the place of the regular wall.
        ExtrusionLoop *top_pass = new ExtrusionLoop(std::move(wall_loop.passes.back()));
        delete *wall_loop.slot;
        *wall_loop.slot = top_pass;
    }
}

} // namespace Slic3r
