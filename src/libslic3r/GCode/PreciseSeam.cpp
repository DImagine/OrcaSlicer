// Precise Seam implementation. Design: docs/HLSD/precise-seam.md
#include "PreciseSeam.hpp"
#include "PreciseSeamInternal.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Model.hpp"
#include "SeamPlacer.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/AABBTreeLines.hpp"
#include "libslic3r/Line.hpp"
#include <algorithm>
#include <atomic>
#include <boost/log/trivial.hpp>
#include <optional>
#include <cstddef>
#include <limits>
#include <string>
#include <cassert>
#include <cmath>
#include <tbb/parallel_for.h>
#include <vector>
#include <utility>

namespace Slic3r {
namespace PreciseSeam {

ModifierRegions prepare_modifier_regions(ExPolygons regions)
{
    ModifierRegions cached;
    cached.reserve(regions.size());
    // Each exterior is measured once; its holes remain in the same cached region.
    for (ExPolygon &region : regions)
        cached.emplace_back(std::move(region));
    return cached;
}

ModifierSlices prepare_modifier_slices(std::vector<ExPolygons> slices)
{
    ModifierSlices cached;
    cached.reserve(slices.size());
    // Keep empty layers so callers can index by the original object layer number.
    for (ExPolygons &layer : slices)
        cached.push_back(prepare_modifier_regions(std::move(layer)));
    return cached;
}

PreparedPerimeter::PreparedPerimeter(const Polygon &perimeter) : polygon(perimeter)
{
    const size_t count = polygon.size();
    if (count < 3)
        return;
    for (size_t edge = 0; edge < count; ++edge)
        if (polygon.points[edge] == polygon.points[(edge + 1) % count])
            return;
    bounds = BoundingBox(polygon.points);
    // Close an open polyline explicitly; prepare its storage once for all modifier queries.
    line.points.reserve(count + 1);
    line.points.insert(line.points.end(), polygon.points.begin(), polygon.points.end());
    line.points.push_back(polygon.points.front());
    valid = true;
}

// Import EnforcedBlockedSeamPoint from SeamPlacerImpl namespace for convenience
using SeamPlacerImpl::EnforcedBlockedSeamPoint;

// Squared on-edge tolerance in coordinate units: Clipper truncates cuts to whole units, under 1 unit
// per axis, so under sqrt(2) from the edge (squared < 2); 2.5 keeps a margin. Independent of SCALING_FACTOR.
static constexpr double MACHINE_PRECISION_SQUARED = 2.5;

// Snapping radius for point insertion, shared by every rule that must match insertion. Deliberately in
// scaled units (1 um by default, 10 um on large beds), well above single-precision candidate steps.
static constexpr coord_t TOLERANCE_LINEAR = 1000;
static constexpr coord_t TOLERANCE_SQUARED = TOLERANCE_LINEAR * TOLERANCE_LINEAR;

namespace detail {

// Use the existing clipping tolerance only to recover rounded coordinates, not to
// bridge gaps between intervals: even a small uncovered interval must remain a gap.
static std::optional<double> parameter_on_edge(const Point &point, const Point &a, const Point &b)
{
    // Original vertices need no floating-point projection.
    if (point == a) return 0.;
    if (point == b) return 1.;
    const Vec2d direction = b.cast<double>() - a.cast<double>();
    const Vec2d offset = point.cast<double>() - a.cast<double>();
    const double squared_length = direction.squaredNorm();
    if (squared_length == 0.)
        return std::nullopt;
    const double parameter = std::clamp(offset.dot(direction) / squared_length, 0., 1.);
    if ((offset - parameter * direction).squaredNorm() > MACHINE_PRECISION_SQUARED)
        return std::nullopt;
    return parameter;
}

static std::optional<ClippedEdgeInterval> interval_on_edge(
    const Point &first, const Point &last, size_t edge, const Polygon &perimeter)
{
    const Point &a = perimeter.points[edge];
    const Point &b = perimeter.points[(edge + 1) % perimeter.size()];
    const auto t0 = parameter_on_edge(first, a, b);
    // Reject this edge before projecting the second point when the first is outside.
    if (!t0)
        return std::nullopt;
    const auto t1 = parameter_on_edge(last, a, b);
    if (!t1 || *t0 == *t1)
        return std::nullopt;
    ClippedEdgeInterval interval{edge, *t0, *t1, first, last};
    if (interval.begin > interval.end) {
        std::swap(interval.begin, interval.end);
        std::swap(interval.first, interval.last);
    }
    // Canonical endpoints make adjacent original edges join at their shared vertex.
    if (interval.begin == 0.) interval.first = a;
    if (interval.end == 1.) interval.last = b;
    return interval;
}

// Interior clipping vertices normally retain the exact source coordinates. Only
// the two cut endpoints need projection; a failed sequence tries the next anchor.
static bool append_exact_fragment(const Polyline &fragment, const Polygon &perimeter,
                                  std::vector<ClippedEdgeInterval> &intervals)
{
    const size_t size = fragment.size();
    const size_t count = perimeter.size();
    assert(size >= 3); // bind_fragment() sends only fragments with interior points here.
    // A longer sequence cannot be one exact pass over this contour.
    if (size > count + 2)
        return false;
    for (size_t anchor = 0; anchor < count; ++anchor) {
        if (fragment.points[1] != perimeter.points[anchor])
            continue;
        // Clipping may return either direction, independent of contour winding.
        for (bool forward : {true, false}) {
            const auto next = [count, forward](size_t index) {
                return forward ? (index + 1) % count : (index + count - 1) % count;
            };
            size_t last = anchor;
            bool matches = true;
            for (size_t i = 2; i + 1 < size; ++i) {
                last = next(last);
                if (fragment.points[i] != perimeter.points[last]) {
                    matches = false;
                    break;
                }
            }
            if (!matches)
                continue;
            const size_t first_edge = forward ? (anchor + count - 1) % count : anchor;
            const size_t last_edge = forward ? last : (last + count - 1) % count;
            const auto first = interval_on_edge(fragment.points.front(), fragment.points[1], first_edge, perimeter);
            const auto end = interval_on_edge(fragment.points[size - 2], fragment.points.back(), last_edge, perimeter);
            if (!first || !end)
                continue;
            // Commit only a complete match, so rejected anchors leave no intervals.
            intervals.push_back(*first);
            size_t vertex = anchor;
            for (size_t i = 2; i + 1 < size; ++i) {
                const size_t adjacent = next(vertex);
                const size_t edge = forward ? vertex : adjacent;
                intervals.push_back({edge, 0., 1., perimeter.points[edge], perimeter.points[(edge + 1) % count]});
                vertex = adjacent;
            }
            intervals.push_back(*end);
            return true;
        }
    }
    return false;
}

// Find the first edge once, then follow the contour without restarting a global search.
bool append_projected_fragment(const Polyline &fragment, const Polygon &perimeter,
                               std::vector<ClippedEdgeInterval> &intervals,
                               FragmentBindingFailure &failure)
{
    const size_t original_size = intervals.size();
    bool forward = true;
    const size_t count = perimeter.size();
    for (size_t i = 1; i < fragment.size(); ++i) {
        if (fragment.points[i - 1] == fragment.points[i])
            continue;
        failure.pair_index = i - 1;
        std::optional<ClippedEdgeInterval> matched;
        const auto is_forward = [&](const ClippedEdgeInterval &candidate) {
            const Vec2d movement = fragment.points[i].cast<double>() - fragment.points[i - 1].cast<double>();
            const Vec2d edge = perimeter.points[(candidate.edge + 1) % count].cast<double>() -
                               perimeter.points[candidate.edge].cast<double>();
            return movement.dot(edge) > 0.;
        };
        if (intervals.size() == original_size) {
            failure.reason = "initial edge not found";
            for (size_t edge = 0; edge < count; ++edge) {
                matched = interval_on_edge(fragment.points[i - 1], fragment.points[i], edge, perimeter);
                if (matched) {
                    // The first suitable edge defines orientation; overlapping visits are unsupported.
                    forward = is_forward(*matched);
                    break;
                }
            }
        } else {
            failure.reason = "non-continuous binding";
            const ClippedEdgeInterval &previous = intervals.back();
            const double end = forward ? previous.end : previous.begin;
            matched = interval_on_edge(fragment.points[i - 1], fragment.points[i], previous.edge, perimeter);
            if (matched && is_forward(*matched) != forward)
                matched.reset();
            // The shared point is the previous pair's last point on the same edge, so continuity holds by
            // construction. Reuse its parameter rather than comparing a recomputed one: the two projections
            // may be compiled differently (e.g. FMA contraction) and differ by an ulp on some platforms.
            if (matched)
                (forward ? matched->begin : matched->end) = end;
            // Crossing to the neighbor is allowed only at their actual shared vertex.
            const size_t vertex = forward ? (previous.edge + 1) % count : previous.edge;
            if (!matched && end == (forward ? 1. : 0.) && fragment.points[i - 1] == perimeter.points[vertex]) {
                const size_t edge = forward ? vertex : (previous.edge + count - 1) % count;
                matched = interval_on_edge(fragment.points[i - 1], fragment.points[i], edge, perimeter);
                if (matched && is_forward(*matched) != forward)
                    matched.reset();
                // The shared point is exactly the neighbor's start (forward) or end (backward) vertex,
                // so parameter_on_edge() returns exactly 0 or 1 for it.
                assert(!matched || (forward ? matched->begin == 0. : matched->end == 1.));
            }
        }
        if (!matched) {
            // Never leave a partial fragment or discard earlier successful fragments.
            intervals.resize(original_size);
            return false;
        }
        intervals.push_back(*matched);
    }
    return intervals.size() > original_size;
}

// Exact binding first, projection as the general path; both leave intervals unchanged on failure.
bool bind_fragment(const Polyline &fragment, const Polygon &perimeter,
                   std::vector<ClippedEdgeInterval> &intervals, FragmentBindingFailure &failure)
{
    if (fragment.size() > 2 && append_exact_fragment(fragment, perimeter, intervals))
        return true;
    return append_projected_fragment(fragment, perimeter, intervals, failure);
}

// Rare-case repair, run only after binding has failed: Clipper may place an end cut a unit or two beside
// a vertex, so it lies on neither edge. Move such an end onto that vertex, taken from the fragment's own
// chain, and let the caller retry with the same strict rules, so a wrong move can only fail again.
// Returns false when nothing was changed.
static bool snap_cuts_to_adjacent_vertices(Polyline &fragment, const Polygon &perimeter)
{
    const size_t count = perimeter.size();
    // Strict comparison matches the vertex snapping in insert_point_into_perimeter.
    const auto close = [](const Point &a, const Point &b) {
        return (a - b).cast<double>().squaredNorm() < double(TOLERANCE_SQUARED);
    };
    const auto snap_target = [&](const Point &cut, const Point &neighbour) -> std::optional<Point> {
        if (cut == neighbour)
            return std::nullopt;
        bool neighbour_is_vertex = false;
        bool ambiguous = false;
        std::optional<Point> target;
        // Scan the whole perimeter before deciding: a cut that is a real vertex must be detected anywhere,
        // and an ambiguity must not pre-empt the neighbour rule below.
        for (size_t j = 0; j < count; ++j) {
            if (perimeter[j] == cut)
                return std::nullopt; // A real vertex, not a rounded cut.
            if (perimeter[j] != neighbour)
                continue;
            neighbour_is_vertex = true;
            for (const Point &candidate : {perimeter[(j + count - 1) % count], perimeter[(j + 1) % count]}) {
                if (!close(cut, candidate))
                    continue;
                if (target && *target != candidate)
                    ambiguous = true;
                target = candidate;
            }
        }
        // A cut between two arbitrary points of one edge is not this rounding case.
        if (!neighbour_is_vertex)
            return std::nullopt;
        // The neighbour itself takes precedence, even when both of its chain vertices are close too
        // (edges shorter than the radius on both sides): the cut is then just a rounded copy of it.
        if (close(cut, neighbour))
            return neighbour;
        // Two different edge-start candidates: keep the failure rather than guess.
        if (ambiguous)
            return std::nullopt;
        return target;
    };
    bool changed = false;
    if (fragment.size() >= 2)
        if (const auto target = snap_target(fragment.points.front(), fragment.points[1])) {
            if (*target == fragment.points[1])
                fragment.points.erase(fragment.points.begin());
            else
                fragment.points.front() = *target;
            changed = true;
        }
    if (fragment.size() >= 2)
        if (const auto target = snap_target(fragment.points.back(), fragment.points[fragment.size() - 2])) {
            if (*target == fragment.points[fragment.size() - 2])
                fragment.points.pop_back();
            else
                fragment.points.back() = *target;
            changed = true;
        }
    return changed;
}

// Location fields shared by the failure and recovery markers, enough to find the layer in a saved project.
static std::string fragment_location(const ExtractionContext &context, size_t fragment_index)
{
    const Layer *layer = context.layer;
    const ModelObject *object = layer && layer->object() ? layer->object()->model_object() : nullptr;
    return " object=" + std::to_string(object ? object->id().id : 0) +
           " modifier=" + std::to_string(context.modifier ? context.modifier->id().id : 0) +
           " layer=" + (layer ? std::to_string(layer->id()) : std::string("unknown")) +
           " z=" + (layer ? std::to_string(layer->slice_z) : std::string("unknown")) +
           " fragment=" + std::to_string(fragment_index);
}

bool append_fragment(const Polyline &fragment, const Polygon &perimeter,
                     std::vector<ClippedEdgeInterval> &intervals,
                     const ExtractionContext &context, size_t fragment_index)
{
    FragmentBindingFailure failure;
    if (bind_fragment(fragment, perimeter, intervals, failure))
        return true;
    // Rare-case fallback: repair the end cuts and retry once. A recovery is logged but is not a failure
    // (no failure count, no user warning).
    const char *outcome = nullptr;
    Polyline cleaned = fragment;
    if (snap_cuts_to_adjacent_vertices(cleaned, perimeter)) {
        FragmentBindingFailure retry_failure;
        // Only a contact shorter than the snapping distance may remain: nothing to bind, nothing lost.
        if (cleaned.size() < 2)
            outcome = "contact";
        else if (bind_fragment(cleaned, perimeter, intervals, retry_failure))
            outcome = "bound";
    }
    // A failed fragment shorter than the snapping radius is a contact: insertion would collapse it to one
    // point anyway, so accept it instead of reporting a failure.
    if (outcome == nullptr && fragment.length() < double(TOLERANCE_LINEAR))
        outcome = "contact";
    if (outcome != nullptr) {
        // Same log budget as failures, counted separately; SeamPlacer::init() reports the total.
        if (!context.quiet && (context.warnings == nullptr ||
            context.warnings->recovered_fragments.fetch_add(1, std::memory_order_relaxed) < failed_fragment_log_limit))
            BOOST_LOG_TRIVIAL(warning) << "[PreciseSeamFragmentRecovered] Rare case resolved after a binding failure"
                << " outcome=" << outcome << fragment_location(context, fragment_index)
                << " original_pair=" << failure.pair_index << " original_reason=" << failure.reason
                << " fragment_points=" << fragment.size() << " cleaned_points=" << cleaned.size()
                << " perimeter_points=" << perimeter.size();
        return true;
    }
    if (context.quiet)
        return false;
    // Reserve a log slot atomically before formatting; callers without shared state log every failure.
    if (context.warnings &&
        context.warnings->failed_fragments.fetch_add(1, std::memory_order_relaxed) >= failed_fragment_log_limit)
        return false;
    // Keep a small marker for investigating a saved project, not a full geometry dump.
    BOOST_LOG_TRIVIAL(warning) << "[PreciseSeamIntersectionFailed] Unable to process intersection"
        << fragment_location(context, fragment_index) << " pair=" << failure.pair_index
        << " reason=" << failure.reason << " fragment_points=" << fragment.size()
        << " perimeter_points=" << perimeter.size();
    return false;
}
} // namespace detail

// Finalize one joined strong segment; only Center needs temporary per-edge lengths.
static void prepare_strong_segment(PerimeterSegment &segment, ModelVolumeType mode, bool compare_lengths)
{
    const bool center = mode == ModelVolumeType::PRECISE_SEAM_CENTER;
    std::vector<double> edge_lengths;
    if (center)
        edge_lengths.reserve(segment.edge_indices.size());
    // A single Left/Right segment needs only its endpoint; Center always needs arc length.
    if (center || compare_lengths) {
        for (size_t i = 0; i < segment.edge_indices.size(); ++i) {
            const double length = (segment.polyline.points[i + 1].cast<double>() -
                                   segment.polyline.points[i].cast<double>()).norm();
            segment.length += length;
            if (center)
                edge_lengths.push_back(length);
        }
    }
    if (mode == ModelVolumeType::PRECISE_SEAM_LEFT) {
        segment.strong_target = StrongSeamTarget{segment.polyline.points.front(), segment.begin.edge_index};
        return;
    }
    segment.strong_target = StrongSeamTarget{segment.polyline.points.back(), segment.end.edge_index};
    if (!center)
        return;

    // The total is now known; reuse local lengths to locate its half without new square roots.
    double remaining = segment.length * 0.5;
    for (size_t i = 0; i < edge_lengths.size(); ++i) {
        const double length = edge_lengths[i];
        if (remaining <= length && length > 0.) {
            const Point &a = segment.polyline.points[i];
            const Point &b = segment.polyline.points[i + 1];
            segment.strong_target = StrongSeamTarget{
                a + ((remaining / length) * (b - a).cast<double>()).cast<coord_t>(), segment.edge_indices[i]};
            return;
        }
        remaining -= length;
    }
}

SegmentExtraction extract_perimeter_segments(const PreparedPerimeter &prepared, const ModifierRegions &modifier,
                                             ModelVolumeType mode, const ExtractionContext &context)
{
    using detail::ClippedEdgeInterval;
    SegmentExtraction result;
    if (!prepared.valid) {
        result.valid = false;
        return result;
    }
    const Polygon &perimeter = prepared.polygon;
    const size_t count = perimeter.size();
    std::vector<const ExPolygon*> nearby_regions;
    for (const ModifierRegion &region : modifier)
        if (!region.polygon.empty() && prepared.bounds.overlap(region.bounds))
            nearby_regions.push_back(&region.polygon);
    if (nearby_regions.empty())
        return result;

    // Clip as an open line, without offsets or contact rules: modifiers should cross the perimeter
    // clearly, and the open path keeps self-touching vertices apart.
    std::vector<ClippedEdgeInterval> intervals;
    Polylines fragments;
    // Clip only accepted regions, with their holes still attached to the exterior.
    for (const ExPolygon *region : nearby_regions) {
        Polylines clipped = intersection_pl(prepared.line, *region);
        for (Polyline &fragment : clipped)
            // Clipper returns a single point where a border only touches the line: nothing to bind.
            if (fragment.size() >= 2)
                fragments.push_back(std::move(fragment));
    }
    for (size_t i = 0; i < fragments.size(); ++i)
        if (!detail::append_fragment(fragments[i], perimeter, intervals, context, i))
            ++result.discarded_fragments;

    std::sort(intervals.begin(), intervals.end(), [](const auto &a, const auto &b) {
        if (a.edge != b.edge) return a.edge < b.edge;
        if (a.begin != b.begin) return a.begin < b.begin;
        return a.end < b.end;
    });
    // Union overlaps on the same occurrence of an edge, never across equal coordinates.
    std::vector<ClippedEdgeInterval> merged;
    for (const ClippedEdgeInterval &interval : intervals) {
        // Fragments touching on an edge meet at the same integer point. Comparing that point as well keeps
        // them joined even if its two projections differ by an ulp on some platforms (e.g. FMA contraction).
        if (!merged.empty() && merged.back().edge == interval.edge &&
            (interval.begin <= merged.back().end || interval.first == merged.back().last)) {
            if (interval.end > merged.back().end) {
                merged.back().end = interval.end;
                merged.back().last = interval.last;
            }
        } else
            merged.push_back(interval);
    }
    const auto position = [count](size_t edge, double parameter) {
        return parameter == 1. ? PerimeterPosition{(edge + 1) % count, 0.}
                               : PerimeterPosition{edge, parameter};
    };
    const auto same_position = [](const PerimeterPosition &a, const PerimeterPosition &b) {
        return a.edge_index == b.edge_index && a.parameter == b.parameter;
    };
    for (const ClippedEdgeInterval &interval : merged) {
        const auto begin = position(interval.edge, interval.begin);
        const auto end = position(interval.edge, interval.end);
        if (result.segments.empty() || !same_position(result.segments.back().end, begin)) {
            PerimeterSegment segment;
            segment.begin = begin;
            segment.polyline.points.push_back(interval.first);
            result.segments.push_back(std::move(segment));
        }
        PerimeterSegment &segment = result.segments.back();
        segment.polyline.points.push_back(interval.last);
        segment.edge_indices.push_back(interval.edge);
        segment.end = end;
    }
    if (result.segments.size() > 1 && same_position(result.segments.back().end, result.segments.front().begin)) {
        // Only the artificial cut at vertex zero can join the last and first intervals.
        PerimeterSegment tail = std::move(result.segments.back());
        result.segments.pop_back();
        PerimeterSegment &head = result.segments.front();
        tail.polyline.points.insert(tail.polyline.points.end(), head.polyline.points.begin() + 1, head.polyline.points.end());
        tail.edge_indices.insert(tail.edge_indices.end(), head.edge_indices.begin(), head.edge_indices.end());
        tail.end = head.end;
        head = std::move(tail);
    }
    // A full loop must cover every original edge, not merely have equal endpoint coordinates.
    result.full_containment = merged.size() == count;
    for (size_t i = 0; result.full_containment && i < count; ++i)
        result.full_containment = merged[i].edge == i && merged[i].begin == 0. && merged[i].end == 1.;
    // A touch on an inclined edge leaves a sub-micron gap. Treat the single segment as full containment
    // when the gap is shorter than the snapping radius, or spans one vertex with both ends within it
    // (exact for edges of 2 um and longer). The 2 um end distance is a cheap pre-filter.
    if (!result.full_containment && result.segments.size() == 1) {
        const PerimeterSegment &only = result.segments.front();
        const Point &first = only.polyline.points.front();
        const Point &last = only.polyline.points.back();
        const double limit = double(TOLERANCE_LINEAR);
        if ((first - last).cast<double>().squaredNorm() < 4. * limit * limit) {
            const auto edge_length = [&perimeter, count](size_t edge) {
                return (perimeter[(edge + 1) % count] - perimeter[edge]).cast<double>().norm();
            };
            // Uncovered length from the segment's end forward to its begin; stops early once it exceeds the limit.
            const auto gap_length = [&](const PerimeterPosition &from, const PerimeterPosition &to) {
                if (from.edge_index == to.edge_index && from.parameter <= to.parameter)
                    return (to.parameter - from.parameter) * edge_length(from.edge_index);
                double gap = (1. - from.parameter) * edge_length(from.edge_index);
                for (size_t edge = (from.edge_index + 1) % count; gap < limit && edge != to.edge_index; edge = (edge + 1) % count)
                    gap += edge_length(edge);
                return gap + to.parameter * edge_length(to.edge_index);
            };
            // Strict comparisons match the vertex snapping in insert_point_into_perimeter.
            const auto within_snap = [limit](const Point &a, const Point &b) {
                return (a - b).cast<double>().squaredNorm() < limit * limit;
            };
            const Point &vertex = perimeter[only.begin.edge_index];
            const bool around_one_vertex = (only.end.edge_index + 1) % count == only.begin.edge_index &&
                                           within_snap(last, vertex) && within_snap(first, vertex);
            if (around_one_vertex || gap_length(only.end, only.begin) < limit)
                result.full_containment = true;
        }
    }
    // Strong skips full containment and weak handles it by type; neither needs strong data for it.
    if (!result.full_containment && is_precise_seam_strong(mode))
        for (PerimeterSegment &segment : result.segments)
            prepare_strong_segment(segment, mode, result.segments.size() > 1);
    return result;
}

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

void init_precise_seam_data(
    std::vector<const ModelVolume*>& strong_volumes_out,
    std::vector<const ModelVolume*>& weak_volumes_out,
    bool& has_strong_out,
    const ModelObject* model_object)
{
    // Clear output vectors
    strong_volumes_out.clear();
    weak_volumes_out.clear();

    if (model_object == nullptr) {
        has_strong_out = false;
        return;
    }

    // Collect and categorize precise seam modifiers
    for (const ModelVolume* volume : model_object->volumes) {
        if (volume->is_precise_seam()) {
            // Strong (Center/Left/Right) and weak (Enforced/Blocked/Neutral) go to separate lists.
            if (volume->is_precise_seam_strong())
                strong_volumes_out.push_back(volume);
            else
                weak_volumes_out.push_back(volume);
        }
    }

    has_strong_out = !strong_volumes_out.empty();

    // Collection already preserves model order, with higher-priority strong modifiers first.
    // Weak modifiers use last-write-wins, so apply the higher-priority ones last.
    std::reverse(weak_volumes_out.begin(), weak_volumes_out.end());
}

// Inserts a point on edge `edge_start_idx`, or reuses a vertex closer than TOLERANCE_LINEAR.
// Returns {final coordinates, index}. Requires at least three vertices (callers ensure it).
static std::pair<Point, size_t> insert_point_into_perimeter(
    const Point &point,
    size_t edge_start_idx,
    Polygon &perimeter_polygon
)
{
    assert(perimeter_polygon.points.size() >= 3);
    size_t perim_max = perimeter_polygon.points.size();

    // Determine edge start and end
    size_t vtx_start = edge_start_idx;
    size_t vtx_end = (edge_start_idx + 1) % perim_max;

    const Point &perim_p_start = perimeter_polygon.points[vtx_start];
    const Point &perim_p_end = perimeter_polygon.points[vtx_end];

    // Check proximity to edge vertices
    coord_t dist_sq_start = (point - perim_p_start).squaredNorm();
    if (dist_sq_start < TOLERANCE_SQUARED) {
        return std::make_pair(perim_p_start, vtx_start);
    }

    coord_t dist_sq_end = (point - perim_p_end).squaredNorm();
    if (dist_sq_end < TOLERANCE_SQUARED) {
        return std::make_pair(perim_p_end, vtx_end);
    }

    // On the closing edge, append instead of inserting at index 0: callers insert in descending source
    // order and rely on earlier indices staying valid.
    size_t insert_pos;
    if (edge_start_idx == perim_max - 1) {
        // Last edge: add to end of vector
        perimeter_polygon.points.push_back(point);
        insert_pos = perimeter_polygon.points.size() - 1;
    } else {
        // Regular edge: insert before end vertex
        insert_pos = vtx_end;
        perimeter_polygon.points.insert(
            perimeter_polygon.points.begin() + insert_pos,
            point
        );
    }

    return std::make_pair(perimeter_polygon.points[insert_pos], insert_pos);
}

// Adds a helper point TOLERANCE_LINEAR after (+1) or before (-1) vertex `point_idx`; skips edges shorter
// than 2 * TOLERANCE_LINEAR. Requires at least three vertices.
static void refine_at_vertex(
    size_t point_idx,
    int direction,
    Polygon &perimeter_polygon
)
{
    assert(perimeter_polygon.points.size() >= 3);
    assert(direction == 1 || direction == -1);

    size_t perim_max = perimeter_polygon.points.size();

    // Determine target edge based on direction
    size_t edge_start_idx, edge_end_idx;

    if (direction == 1) {
        // Direction +1: insertion AFTER point_idx (edge point_idx → point_idx+1)
        edge_start_idx = point_idx;
        edge_end_idx = (point_idx + 1) % perim_max;
    } else {
        // Direction -1: insertion BEFORE point_idx (edge point_idx-1 → point_idx)
        edge_start_idx = (point_idx + perim_max - 1) % perim_max;
        edge_end_idx = point_idx;
    }

    const Point &edge_start = perimeter_polygon.points[edge_start_idx];
    const Point &edge_end = perimeter_polygon.points[edge_end_idx];

    // Calculate edge length
    Vec2d edge_vector = (edge_end - edge_start).cast<double>();
    double edge_length = edge_vector.norm();

    // The helper needs TOLERANCE_LINEAR of edge on both sides.
    if (edge_length < 2.0 * TOLERANCE_LINEAR) {
        return;  // Edge too short - new point would be too close to end
    }

    // Calculate new point coordinates: edge_start + TOLERANCE_LINEAR * direction_normalized
    Vec2d direction_normalized = edge_vector / edge_length;
    // Place helper point near point_idx: after it for +1, before it for -1
    auto offset = (TOLERANCE_LINEAR * direction_normalized).cast<coord_t>();
    Point new_point = (direction == 1)
        ? Point(edge_start + offset)
        : Point(edge_end   - offset);

    // On the closing edge, append instead of inserting at index 0: callers insert in descending source
    // order and rely on earlier indices staying valid.
    if (edge_start_idx == perim_max - 1) {
        // Last edge: add to end of vector
        perimeter_polygon.points.push_back(new_point);
    } else {
        // Regular edge: insert before end vertex
        perimeter_polygon.points.insert(
            perimeter_polygon.points.begin() + edge_end_idx,
            new_point
        );
    }
}

// The type to apply on this layer, or nullopt where the rules exclude the layer. Without rules
// (`effective` is null) every modifier keeps its own type on every layer.
static std::optional<ModelVolumeType> active_type(const EffectiveModifiers *effective, const ModelVolume *modifier, size_t layer_id)
{
    if (effective != nullptr)
        if (const auto it = effective->find(modifier); it != effective->end())
            return it->second.active_on(layer_id) ? std::optional(it->second.type) : std::nullopt;
    return modifier->type();
}

// The magnet of a modifier with layer rules, or null.
static const EffectiveModifier *magnet_of(const EffectiveModifiers *effective, const ModelVolume *modifier)
{
    if (effective != nullptr)
        if (const auto it = effective->find(modifier); it != effective->end() && it->second.magnet)
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

// Strong priority is per modifier, never a global maximum across different modifiers.
std::optional<Point> insert_strong_seam_point(
    const std::vector<const ModelVolume*> &strong_volumes,
    Polygon &polygon,
    const PreparedPerimeter &prepared,
    const Layer *layer,
    const ModifierRegionsCache &slices_cache,
    PreciseSeamWarnings* warnings,
    const EffectiveModifiers *effective,
    bool quiet)
{
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

// Convert ModelVolumeType of weak modifier to EnforcedBlockedSeamPoint.
// Precondition: a weak precise-seam type; init_precise_seam_data() puts only those in the weak list.
// Exhaustive switch (no default) so -Wswitch flags any future PRECISE_SEAM_* additions.
static EnforcedBlockedSeamPoint convert_weak_modifier_type(ModelVolumeType type) {
    assert(is_precise_seam_weak(type));
    switch (type) {
        case ModelVolumeType::PRECISE_SEAM_ENFORCED:
            return EnforcedBlockedSeamPoint::Enforced;
        case ModelVolumeType::PRECISE_SEAM_BLOCKED:
            return EnforcedBlockedSeamPoint::Blocked;
        case ModelVolumeType::PRECISE_SEAM_NEUTRAL:
            return EnforcedBlockedSeamPoint::Neutral;
        // Non-weak types are unreachable by precondition; listed to keep the switch exhaustive.
        case ModelVolumeType::INVALID:
        case ModelVolumeType::MODEL_PART:
        case ModelVolumeType::NEGATIVE_VOLUME:
        case ModelVolumeType::PARAMETER_MODIFIER:
        case ModelVolumeType::SUPPORT_BLOCKER:
        case ModelVolumeType::SUPPORT_ENFORCER:
        case ModelVolumeType::PRECISE_SEAM_CENTER:
        case ModelVolumeType::PRECISE_SEAM_LEFT:
        case ModelVolumeType::PRECISE_SEAM_RIGHT:
            break;
    }
    assert(false && "convert_weak_modifier_type called with non-weak type");
    return EnforcedBlockedSeamPoint::Neutral;
}

// Consume ready boundaries in modifier priority order; geometry extraction is separate.
static std::vector<WeakModifierSegment> prepare_weak_modifier_segments(
    std::vector<WeakModifierSegment> result, Polygon &polygon)
{
    // If no segments, return empty vector
    if (result.empty()) {
        return result;
    }

    // Source positions have the same order as arc length, without measuring the perimeter.
    struct PointToInsert {
        size_t segment_idx;
        bool is_left;
        PerimeterPosition position;
    };
    std::vector<PointToInsert> points_to_insert;
    points_to_insert.reserve(result.size() * 2);
    for (size_t seg_idx = 0; seg_idx < result.size(); ++seg_idx) {
        const WeakModifierSegment &seg = result[seg_idx];
        // A whole-perimeter zone has no boundaries to insert.
        if (seg.whole_perimeter)
            continue;
        points_to_insert.push_back({seg_idx, true, seg.left_position});
        points_to_insert.push_back({seg_idx, false, seg.right_position});
    }

    // Descending order preserves pending source indices; vertex zero is canonical (0, 0).
    // Only insertion events are reordered: modifier priority in result remains unchanged.
    std::sort(points_to_insert.begin(), points_to_insert.end(),
              [](const PointToInsert &a, const PointToInsert &b) {
                  if (a.position.edge_index != b.position.edge_index)
                      return a.position.edge_index > b.position.edge_index;
                  return a.position.parameter > b.position.parameter;
              });

    for (const PointToInsert &pt : points_to_insert) {
        WeakModifierSegment &seg = result[pt.segment_idx];
        Point &point_coords = pt.is_left ? seg.left_point : seg.right_point;
        // Take the final coordinates: a boundary close to an existing vertex snaps to it.
        point_coords = insert_point_into_perimeter(point_coords, pt.position.edge_index, polygon).first;
    }

    // Add helpers at boundary vertices, scanning original vertices backwards so insertions do not shift
    // pending indices.
    for (size_t poly_idx = polygon.size(); poly_idx-- > 0; ) {
        bool refine_before = false;
        bool refine_after = false;
        for (const WeakModifierSegment &segment : result) {
            // Its placeholder points must not match a real vertex.
            if (segment.whole_perimeter)
                continue;
            refine_before |= polygon[poly_idx] == segment.left_point;
            refine_after |= polygon[poly_idx] == segment.right_point;
        }
        // Insert after first: inserting before would shift the current vertex index.
        if (refine_after)
            refine_at_vertex(poly_idx, +1, polygon);
        if (refine_before)
            refine_at_vertex(poly_idx, -1, polygon);
    }

    // Split edges in enforced zones into segments ≤ enforcer_oversampling_distance
    // Determine type pattern for each polygon edge (sequential application of hierarchy)
    std::vector<EnforcedBlockedSeamPoint> edge_types(polygon.size(), EnforcedBlockedSeamPoint::Neutral);

    // Every boundary was inserted or snapped above and later steps only add vertices, so it is found.
    // A linear scan is fine: there are only a few zones per perimeter.
    auto find_point_index = [&](const Point &pt) -> size_t {
        const auto it = std::find(polygon.points.begin(), polygon.points.end(), pt);
        assert(it != polygon.points.end());
        return size_t(it - polygon.points.begin());
    };

    // Last write wins: zones are sorted lowest priority first.
    for (const auto &segment : result) {
        // A whole-perimeter zone types every edge, so Enforced subdivides the whole perimeter, as
        // painting it green all round would.
        if (segment.whole_perimeter) {
            std::fill(edge_types.begin(), edge_types.end(), segment.type);
            continue;
        }
        const size_t left_idx = find_point_index(segment.left_point);
        const size_t right_idx = find_point_index(segment.right_point);
        // Practically unreachable, but a miss would write past edge_types in Release builds: skip the zone.
        if (left_idx == polygon.size() || right_idx == polygon.size()) {
            BOOST_LOG_TRIVIAL(error) << "PreciseSeam: weak boundary not found in perimeter, skipping zone";
            continue;
        }

        // A zone covers edges [left_idx, right_idx): the edge from the right boundary is outside and must
        // not be subdivided. Edge types only drive subdivision; candidates are typed per point later.
        for (size_t idx = left_idx; idx != right_idx; idx = (idx + 1) % polygon.size())
            edge_types[idx] = segment.type;
    }

    // Boundaries and helpers are already in place; without enforced edges nothing is subdivided,
    // so skip rebuilding the polygon.
    if (std::find(edge_types.begin(), edge_types.end(), EnforcedBlockedSeamPoint::Enforced) == edge_types.end())
        return result;

    // Split enforced edges into small segments
    const double STEP_SCALED = scale_(SeamPlacer::enforcer_oversampling_distance);

    // Collect new list of polygon points with split enforced edges
    Points new_points;
    new_points.reserve(polygon.size() * 10); // Approximate estimate

    for (size_t i = 0; i < polygon.size(); ++i) {
        size_t next_i = (i + 1) % polygon.size();
        Point p_start = polygon[i];
        Point p_end = polygon[next_i];

        // Add current vertex
        new_points.push_back(p_start);

        // Check edge type [i, next_i]
        if (edge_types[i] != EnforcedBlockedSeamPoint::Enforced) {
            continue; // Not enforced - don't split
        }

        // Calculate subdivision parameters
        Vec2d edge_vec = (p_end - p_start).cast<double>();
        double edge_length = edge_vec.norm();

        if (edge_length <= STEP_SCALED) {
            continue; // Edge too short - don't split
        }

        // Number of segments: ceil(L / S)
        size_t num_segments = static_cast<size_t>(std::ceil(edge_length / STEP_SCALED));

        // Uniform step: L / num_segments (each segment ≤ S)
        double actual_step = edge_length / num_segments;
        Vec2d step_vec = edge_vec.normalized() * actual_step;

        // Add intermediate points incrementally
        Vec2d current_pos = p_start.cast<double>();
        for (size_t j = 1; j < num_segments; ++j) {
            current_pos += step_vec;
            new_points.push_back(current_pos.cast<coord_t>());
        }
    }

    // Replace polygon points with new ones (with split enforced edges)
    polygon.points = std::move(new_points);

    return result;
}

// Collect every modifier's ready segments before insertions can shift source edge indices.
std::vector<WeakModifierSegment> collect_weak_modifier_segments(
    const std::vector<const ModelVolume*> &weak_volumes,
    Polygon &polygon,
    const PreparedPerimeter &prepared,
    const Layer *layer,
    const ModifierRegionsCache &slices_cache,
    PreciseSeamWarnings* warnings,
    const EffectiveModifiers *effective)
{
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

void apply_weak_modifiers_to_perimeter(
    const std::vector<WeakModifierSegment> &weak_segments,
    PrintObjectSeamData::LayerSeams &result,
    const SeamPlacerImpl::Perimeter &perimeter,
    bool &some_point_enforced,
    PreciseSeamWarnings *warnings)
{
    // Get z-coordinate for unscaling boundary points
    const float z_coord = result.points[perimeter.start_index].position.z();
    const size_t perimeter_size = perimeter.end_index - perimeter.start_index;

    // Helper lambda: search for point index in result.points by unscaled coordinates.
    // Linear scan is intentional — O(N×M) is acceptable for typical M ≤ 5 weak segments.
    auto find_point_index = [&](const Point &pt) -> std::optional<size_t> {
        Vec2f unscaled_pt = unscale(pt).cast<float>();
        Vec3f target(unscaled_pt.x(), unscaled_pt.y(), z_coord);
        for (size_t i = perimeter.start_index; i < perimeter.end_index; ++i) {
            if (result.points[i].position == target) return i - perimeter.start_index;
        }
        return std::nullopt;
    };

    // Apply weak modifiers sequentially: sorted low-priority-first,
    // so higher-priority modifiers (higher in object tree) overwrite via last-write-wins.
    for (size_t seg_idx = 0; seg_idx < weak_segments.size(); ++seg_idx) {
        const auto &segment = weak_segments[seg_idx];

        // A whole-perimeter zone types every candidate, including painting's oversampled points.
        if (segment.whole_perimeter) {
            for (size_t i = perimeter.start_index; i < perimeter.end_index; ++i)
                result.points[i].type = segment.type;
            if (segment.type == EnforcedBlockedSeamPoint::Enforced)
                some_point_enforced = true;
            continue;
        }

        // Find boundary point indices in result.points
        std::optional<size_t> left_idx = find_point_index(segment.left_point);
        std::optional<size_t> right_idx = find_point_index(segment.right_point);

        if (!left_idx.has_value() || !right_idx.has_value()) {
            BOOST_LOG_TRIVIAL(error) << "PreciseSeam: boundary point not found in perimeter, skipping segment";
            continue;
        }

        // Apply modifier type to range [left_idx, right_idx] with wraparound
        for (size_t idx = left_idx.value(); ; idx = (idx + 1) % perimeter_size) {
            result.points[perimeter.start_index + idx].type = segment.type;
            if (idx == right_idx.value()) break;
        }

        if (segment.type == EnforcedBlockedSeamPoint::Enforced) {
            some_point_enforced = true;
        }
    }

    // Only the highest whole zone is in effect, and a higher band may still unblock part of it.
    if (warnings == nullptr)
        return;
    const auto last_whole = std::find_if(weak_segments.rbegin(), weak_segments.rend(),
                                         [](const WeakModifierSegment &segment) { return segment.whole_perimeter; });
    if (last_whole == weak_segments.rend() || last_whole->type != EnforcedBlockedSeamPoint::Blocked)
        return;
    for (size_t i = perimeter.start_index; i < perimeter.end_index; ++i)
        if (result.points[i].type != EnforcedBlockedSeamPoint::Blocked)
            return;
    if (const auto usage = warnings->modifier_usage.find(last_whole->modifier); usage != warnings->modifier_usage.end())
        PreciseSeamWarnings::set(usage->second.fully_blocked);
}

// Restore precise seam positions that may have been modified
void restore_precise_seam_positions(std::vector<PrintObjectSeamData::LayerSeams> &layers) {
  using SeamPlacerImpl::SeamCandidate;
  using SeamPlacerImpl::Perimeter;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, layers.size()),
    [&layers](tbb::blocked_range<size_t> r) {
      for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
        std::vector<SeamCandidate> &layer_perimeter_points = layers[layer_idx].points;
        // Iterate over perimeters (jump by end_index)
        for (size_t current = 0; current < layer_perimeter_points.size();
             current = layer_perimeter_points[current].perimeter.end_index) {
          Perimeter &perimeter = layer_perimeter_points[current].perimeter;
          if (perimeter.precise_seam_point.has_value()) {
            perimeter.final_seam_position = perimeter.precise_seam_point.value();
            perimeter.seam_index = perimeter.precise_seam_index;
          }
        }
      }
    });
}

} // namespace PreciseSeam
} // namespace Slic3r
