#include "OuterWallPasses.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
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

// Point of the external loop the passes are derived from.
struct LoopSample
{
    Point  point;
    // Unit vector pointing out of the material, scaled so that moving the point along it by d offsets the adjacent
    // edges by d (miter offset of a polygon vertex).
    Vec2d  outward;
    // Offset of each pass outline from the middle outline at this point, positive outwards, scaled.
    std::vector<double> offsets;
    // Whether this is a vertex of the loop or a point added inside an edge.
    bool   vertex;
};

// Build the passes of one external loop. Returns false if the loop has to be printed the regular way.
bool make_passes(const ExtrusionLoop                 &loop,
                 const OuterWallPassSlices           &slices,
                 const OutlineDistancer              &mid_outline,
                 const std::vector<OutlineDistancer> &pass_outlines,
                 const Flow                          &pass_flow,
                 bool                                 with_holes,
                 std::vector<ExtrusionLoop>          &passes_out)
{
    const double nominal_width = slices.nominal_width();
    // Only a plain outer wall of the nominal width is split: overhanging parts keep their bridging flow and
    // a wall held back until the infill is printed must not be printed first.
    if (loop.print_after_infill)
        return false;
    for (const ExtrusionPath &path : loop.paths)
        if (path.role() != erExternalPerimeter || std::abs(double(path.width) - nominal_width) > 1e-3)
            return false;

    Polygon polygon = loop.polygon();
    polygon.remove_duplicate_points();
    const Points &points = polygon.points;
    if (points.size() < 3)
        return false;

    // Sample the loop densely enough to follow a slope changing along a straight edge.
    const double            max_segment = scaled<double>(0.5);
    const size_t            num_passes  = slices.count();
    std::vector<LoopSample> samples;
    std::vector<Vec2d>      edge_normals;
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
            samples.emplace_back(std::move(sample));
        }
    }

    // Measure the offsets of the pass outlines and find out, which side of the loop the material is on.
    const double min_offset = scaled<double>(slices.min_width - nominal_width);
    const double max_offset = scaled<double>(slices.max_width - nominal_width);
    // The material is inside a contour and outside a hole. The normals point to the right of the loop.
    const bool is_hole = (loop.loop_role() & elrHole) == elrHole;
    if (! with_holes && is_hole)
        return false;
    const double side = polygon.is_counter_clockwise() != is_hole ? 1. : -1.;
    // The offsets are measured against the nearest edge of the outlines, which has to be the edge the loop follows:
    // the distance from the middle outline, positive inside, has to fall towards the outer side of the loop.
    // It is probed on both sides of the loop, which works even for a loop running along the outline itself, as it
    // does with an XY size compensation of half the wall width. A loop the compensations moved beside another edge
    // of the outlines, such as the hole of a thin wall they thickened, is left as it is.
    const double side_probe     = scaled<double>(0.05);
    bool         side_confirmed = false;
    for (LoopSample &sample : samples) {
        double mid_distance;
        Vec2d  mid_nearest;
        if (! signed_distance(mid_outline, sample.point, mid_distance, mid_nearest))
            return false;
        const Vec2d right = side_probe * sample.outward.normalized();
        double      distance_right;
        double      distance_left;
        if (signed_distance(mid_outline, Point((sample.point.cast<double>() + right).cast<coord_t>()), distance_right, mid_nearest) &&
            signed_distance(mid_outline, Point((sample.point.cast<double>() - right).cast<coord_t>()), distance_left, mid_nearest) &&
            std::abs(distance_left - distance_right) > 0.5 * side_probe) {
            if ((distance_left > distance_right) != (side > 0.))
                return false;
            side_confirmed = true;
        }
        sample.offsets.reserve(num_passes);
        for (const OutlineDistancer &pass_outline : pass_outlines) {
            double pass_distance;
            Vec2d  pass_nearest;
            if (! signed_distance(pass_outline, sample.point, pass_distance, pass_nearest))
                return false;
            sample.offsets.emplace_back(std::clamp(pass_distance - mid_distance, min_offset, max_offset));
        }
    }
    if (! side_confirmed)
        return false;

    // thick_polyline_to_multi_path() takes the width of the rectangle of the same cross section as the extrusion.
    const double width_correction = scaled<double>(pass_flow.height() * (1. - 0.25 * PI));
    const double linear_tolerance = scaled<double>(0.005);
    const float  width_tolerance  = scaled<float>(0.02);
    passes_out.clear();
    passes_out.reserve(num_passes);
    for (size_t pass = 0; pass < num_passes; ++ pass) {
        ThickPolyline polyline;
        std::vector<coordf_t> widths;
        auto append = [&](const LoopSample &sample) {
            const double offset = sample.offsets[pass];
            polyline.points.emplace_back((sample.point.cast<double>() + side * 0.5 * offset * sample.outward).cast<coord_t>());
            widths.emplace_back(scaled<double>(nominal_width) + offset - width_correction);
        };
        for (size_t i = 0; i < samples.size();) {
            // Samples of one edge: [i, edge_end), the next vertex is at edge_end.
            size_t edge_end = i + 1;
            while (edge_end < samples.size() && ! samples[edge_end].vertex)
                ++ edge_end;
            append(samples[i]);
            // Keep the points inside the edge only if the offset does not change linearly along it.
            const double offset_a = samples[i].offsets[pass];
            const double offset_b = samples[edge_end % samples.size()].offsets[pass];
            bool         linear   = true;
            for (size_t j = i + 1; j < edge_end && linear; ++ j) {
                const double t = double(j - i) / double(edge_end - i);
                linear = std::abs(samples[j].offsets[pass] - (offset_a + t * (offset_b - offset_a))) <= linear_tolerance;
            }
            if (! linear)
                for (size_t j = i + 1; j < edge_end; ++ j)
                    append(samples[j]);
            i = edge_end;
        }
        // Close the loop.
        polyline.points.emplace_back(polyline.points.front());
        widths.emplace_back(widths.front());
        polyline.width.reserve(2 * (polyline.points.size() - 1));
        for (size_t i = 0; i + 1 < polyline.points.size(); ++ i) {
            polyline.width.emplace_back(widths[i]);
            polyline.width.emplace_back(widths[i + 1]);
        }

        // Split the pass into paths of a constant width. The steps are finer than those of thin walls
        // to keep the outer surface smooth.
        ExtrusionPaths paths = thick_polyline_to_multi_path(polyline, erExternalPerimeter, pass_flow, width_tolerance, width_tolerance).paths;
        if (paths.empty() || paths.front().first_point() != paths.back().last_point())
            return false;
        ExtrusionLoop pass_loop(std::move(paths), loop.loop_role());
        pass_loop.inset_idx       = loop.inset_idx;
        pass_loop.outer_wall_pass = true;
        passes_out.emplace_back(std::move(pass_loop));
    }
    return true;
}

void split_loops_recursive(ExtrusionEntityCollection              &collection,
                           std::vector<ExtrusionEntityCollection> &lower_passes,
                           const OuterWallPassSlices              &slices,
                           const OutlineDistancer                 &mid_outline,
                           const std::vector<OutlineDistancer>    &pass_outlines,
                           const Flow                             &pass_flow,
                           bool                                    with_holes)
{
    std::vector<ExtrusionLoop> passes;
    for (ExtrusionEntity *&entity : collection.entities) {
        if (entity->is_collection()) {
            split_loops_recursive(*static_cast<ExtrusionEntityCollection*>(entity), lower_passes, slices, mid_outline, pass_outlines,
                                  pass_flow, with_holes);
        } else if (entity->is_loop() && entity->inset_idx == 0 &&
                   make_passes(*static_cast<const ExtrusionLoop*>(entity), slices, mid_outline, pass_outlines, pass_flow, with_holes, passes)) {
            for (size_t pass = 0; pass + 1 < passes.size(); ++ pass)
                lower_passes[pass].entities.emplace_back(new ExtrusionLoop(std::move(passes[pass])));
            // The top pass is printed at the layer height and takes the place of the regular outer wall.
            delete entity;
            entity = new ExtrusionLoop(std::move(passes.back()));
        }
    }
}

} // namespace

bool outer_wall_passes_fit(const ExPolygons              &mid,
                           const std::vector<ExPolygons> &passes,
                           double                         min_width,
                           double                         max_width,
                           bool                           with_holes)
{
    const double nominal_width = 0.5 * (min_width + max_width);
    // Slices taken at different heights never match exactly: the contours are simplified with this precision.
    const float  tolerance     = scaled<float>(0.003);
    // Slivers the polygon clipping leaves along edges that coincide. Kept far below the tolerance, as whatever
    // the opening removes adds to it.
    const float  sliver        = scaled<float>(0.0005);
    const ExPolygons mid_outline = comparable_outline(mid, with_holes);
    // The widest pass may reach this far out and the narrowest pass has to reach at least this far.
    const Polygons   outer_limit = offset(mid_outline, scaled<float>(max_width - nominal_width) + tolerance);
    const ExPolygons inner_limit = offset_ex(mid_outline, - scaled<float>(nominal_width - min_width) - tolerance);
    for (const ExPolygons &pass : passes) {
        const ExPolygons pass_outline = comparable_outline(pass, with_holes);
        if (! opening_ex(diff_ex(pass_outline, outer_limit), sliver).empty() ||
            ! opening_ex(diff_ex(inner_limit, pass_outline), sliver).empty())
            return false;
    }
    return true;
}

void split_outer_wall_loops(ExtrusionEntityCollection              &perimeters,
                            std::vector<ExtrusionEntityCollection> &lower_passes,
                            const OuterWallPassSlices              &slices,
                            const Flow                             &pass_flow,
                            bool                                    with_holes)
{
    if (slices.count() < 2)
        return;
    lower_passes.resize(slices.count() - 1);
    const OutlineDistancer        mid_outline(to_lines(slices.mid));
    std::vector<OutlineDistancer> pass_outlines;
    pass_outlines.reserve(slices.count());
    for (const ExPolygons &pass : slices.passes)
        pass_outlines.emplace_back(to_lines(pass));
    split_loops_recursive(perimeters, lower_passes, slices, mid_outline, pass_outlines, pass_flow, with_holes);
}

} // namespace Slic3r
