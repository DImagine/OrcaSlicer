#include "PreciseSeam.hpp"
#include "PreciseSeamInternal.hpp"
#include "SeamPlacer.hpp"
#include "libslic3r/BoundingBox.hpp"
#include <algorithm>
#include <boost/log/trivial.hpp>
#include <tbb/parallel_for.h>

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

// Machine precision for checking exact coordinate matching (squared distance)
// Ideally, intersection points should match perimeter vertices bitwise,
// but we account for possible machine rounding errors in Clipper calculations
// Actual deviations: maximum ~0.27, using 2.5 with margin (nanometers)
static constexpr double MACHINE_PRECISION_SQUARED = 2.5;

// Tolerance for checking proximity when inserting seam points into perimeter
static const coord_t TOLERANCE_LINEAR = scale_(0.001);  // 1.0 micrometers
static const coord_t TOLERANCE_SQUARED = TOLERANCE_LINEAR * TOLERANCE_LINEAR;

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
    if (size < 3 || size > count + 2)
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
                if (matched && (is_forward(*matched) != forward ||
                                (forward ? matched->begin != 0. : matched->end != 1.)))
                    matched.reset();
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

// Safety net for a rare case: the modifier boundary must cross the perimeter within about a
// nanometre of a source vertex (Clipper is deterministic, so a prismatic model may repeat it on
// every layer). Not part of the normal binding path; used only after binding has failed.
// Clipper may place a cut at a vertex's height but a few nanometres to the side of it
// (scanbeam clamping, X taken from the modifier edge). The end pair then either collapses to the
// vertex parameter (rejected as zero-length) or misses both neighbouring edges, and the whole
// fragment is lost. Reproduced on OrcaSlicer's own Clipper by a randomized search with borders
// passing within 5 nm of a vertex: about 1 fragment in 20 000 of those needed this cleanup.
//
// The end cut is replaced by a vertex of the fragment's own chain, never by a nearest vertex found
// elsewhere on the perimeter. Candidates are only the cut's neighbour in the fragment, when that is
// a source vertex (the cut is a rounded copy of it and is dropped), or a vertex sharing a source
// edge with that neighbour (the cut stands for the start of that edge). Two different candidates
// (an edge shorter than the radius, or a self-touching contour) leave the fragment unchanged, and so
// does an end that is itself a source vertex (the open line's start/end or an exact cut): it is not
// a rounded cut, and moving it could drop a real edge shorter than the radius.
// The radius is TOLERANCE_LINEAR, the distance at which boundary insertion snaps to an existing
// vertex anyway, so the result matches a successful binding. The caller re-binds with the usual
// strict direction and continuity rules, so a wrong candidate can only fail, never bind elsewhere.
// Returns false when nothing was changed, so the caller does not retry.
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
        std::optional<Point> target;
        // Every visit of the neighbour's coordinate contributes its two chain-adjacent vertices.
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
                    return std::nullopt; // Ambiguous: keep the failure rather than guess.
                target = candidate;
            }
        }
        // A cut between two arbitrary points of one edge is not this rounding case.
        if (!neighbour_is_vertex)
            return std::nullopt;
        // The neighbour itself takes precedence: the cut is then just a rounded copy of it.
        if (close(cut, neighbour))
            return neighbour;
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
    // Fallback for a rare rounding case: retry the same binding rules once on the cleaned fragment.
    // The normal path never gets here. A recovery is not a failure (no failure count, no user warning),
    // but it leaves its own marker so that any later problem can be traced to it.
    // A full failure falls through to the usual failure marker with the original reason.
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
    // Any failed fragment shorter than the snapping distance is a contact as well, for example two
    // cuts in the middle of an edge that project to one parameter: insertion would collapse it onto
    // one point anyway. Accept it as clipping returned it instead of reporting a failure.
    if (outcome == nullptr && fragment.length() < double(TOLERANCE_LINEAR))
        outcome = "contact";
    if (outcome != nullptr) {
        // Same log budget as failures, counted separately; SeamPlacer::init() reports the total.
        if (context.warnings == nullptr ||
            context.warnings->recovered_fragments.fetch_add(1, std::memory_order_relaxed) < failed_fragment_log_limit)
            BOOST_LOG_TRIVIAL(warning) << "[PreciseSeamFragmentRecovered] Rare case resolved after a binding failure"
                << " outcome=" << outcome << fragment_location(context, fragment_index)
                << " original_pair=" << failure.pair_index << " original_reason=" << failure.reason
                << " fragment_points=" << fragment.size() << " cleaned_points=" << cleaned.size()
                << " perimeter_points=" << perimeter.size();
        return true;
    }
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

    // Accept the clipper's boundary behavior without offsets or special contact handling;
    // modifiers should cross the perimeter unambiguously.
    // The open overload avoids coordinate-only recombination at self-touching vertices.
    std::vector<ClippedEdgeInterval> intervals;
    Polylines fragments;
    // Clip only accepted regions, with their holes still attached to the exterior.
    for (const ExPolygon *region : nearby_regions) {
        Polylines clipped = intersection_pl(prepared.line, *region);
        for (Polyline &fragment : clipped)
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
    // A boundary that only touches the perimeter is full containment by policy. On an axis-aligned
    // edge or at a vertex clipping splits the line exactly at the touch and the check above holds. On
    // an inclined edge the touching point is usually not representable on the integer grid, so the
    // boundary pokes a few nm across and leaves a real gap; the single segment then covers everything
    // but that gap. Weak insertion would snap both its boundaries onto one vertex (1 um radius) and turn
    // the intended zone into one candidate, and strong would put the seam at the touch. Treat such a
    // gap as full containment instead, exactly when insertion would collapse it:
    // - the uncovered length is below the snapping distance (a gap inside one edge collapses onto the
    //   boundary inserted first), or
    // - the gap spans exactly one vertex and both ends lie within the snapping distance of it (each
    //   end snaps onto that vertex from its own edge, even if the gap itself is up to 2 um long).
    // Stage 1 is a cheap filter on the normal path: both cases bring the ends closer than 2 um. Narrow
    // bands, outside contacts, sharp spikes and nearly touching contour parts can pass it, but their
    // uncovered part is long and not around one vertex, so they stay ordinary segments (a band below
    // 1 um then becomes a single weak candidate, which is its expected result).
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
    // Full containment is skipped by both consumers; retain geometry but prepare no strong data.
    if (!result.full_containment && is_precise_seam_strong(mode))
        for (PerimeterSegment &segment : result.segments)
            prepare_strong_segment(segment, mode, result.segments.size() > 1);
    return result;
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
            ModelVolumeType type = volume->type();
            // Categorization: strong modifiers have priority
            if (type == ModelVolumeType::PRECISE_SEAM_CENTER ||
                type == ModelVolumeType::PRECISE_SEAM_LEFT ||
                type == ModelVolumeType::PRECISE_SEAM_RIGHT) {
                strong_volumes_out.push_back(volume);
            } else {
                // ENFORCED, BLOCKED, NEUTRAL - weak modifiers (processed later)
                weak_volumes_out.push_back(volume);
            }
        }
    }

    has_strong_out = !strong_volumes_out.empty();

    // Collection already preserves model order, with higher-priority strong modifiers first.
    // Weak modifiers use last-write-wins, so apply the higher-priority ones last.
    std::reverse(weak_volumes_out.begin(), weak_volumes_out.end());
}

// Insert point into perimeter with proximity check to existing vertices
// If point is close to vertex (< TOLERANCE_SQUARED) - use existing vertex
// Returns pair: {final coordinates, point index in polygon}
// edge_start_idx is start vertex of edge containing point
static std::optional<std::pair<Point, size_t>> insert_point_into_perimeter(
    const Point &point,
    size_t edge_start_idx,
    Polygon &perimeter_polygon
)
{
    // Check input data
    if (perimeter_polygon.points.size() < 3) {
        return std::nullopt;  // Polygon must be at least a triangle
    }

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

    // Insert point into perimeter
    // IMPORTANT: Special handling for the last edge to preserve indexing for subsequent insertions.
    // If this is the last edge (edge_start_idx == perim_max - 1), we append to the end instead of
    // inserting at position 0 (which would shift all indices). This allows sorting points by
    // descending source position and inserting them without invalidating previously computed indices.
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

// Insert new point at distance TOLERANCE_LINEAR from specified perimeter vertex
// Insertion direction specified by direction parameter: +1 = after vertex, -1 = before vertex
// If target edge length < 2*TOLERANCE_LINEAR, insertion not performed (new point would be too close to edge end)
// Returns true if point was inserted, false otherwise
// point_idx is index of perimeter vertex from which insertion is performed
static bool refine_at_vertex(
    size_t point_idx,
    int direction,
    Polygon &perimeter_polygon
)
{
    // Check input data
    if (perimeter_polygon.points.size() < 3) {
        return false;  // Polygon must be at least a triangle
    }

    if (direction != 1 && direction != -1) {
        return false;  // Direction must be +1 or -1
    }

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

    // Check if edge is long enough for insertion
    // New point must be at distance TOLERANCE_LINEAR from start
    // and at distance >= TOLERANCE_LINEAR from end
    if (edge_length < 2.0 * TOLERANCE_LINEAR) {
        return false;  // Edge too short - new point would be too close to end
    }

    // Calculate new point coordinates: edge_start + TOLERANCE_LINEAR * direction_normalized
    Vec2d direction_normalized = edge_vector / edge_length;
    // Place helper point near point_idx: after it for +1, before it for -1
    auto offset = (TOLERANCE_LINEAR * direction_normalized).cast<coord_t>();
    Point new_point = (direction == 1)
        ? Point(edge_start + offset)
        : Point(edge_end   - offset);

    // Insert point into perimeter
    // IMPORTANT: Special handling of last edge to preserve indexing for subsequent insertions.
    // If this is last edge (edge_start_idx == perim_max - 1), add point to end of vector
    // instead of inserting at position 0 (which would shift all indices). This allows sorting points
    // by descending source position and inserting them without invalidating previously computed indices.
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

    return true;
}

// Strong priority is per modifier, never a global maximum across different modifiers.
std::optional<Point> insert_strong_seam_point(
    const std::vector<const ModelVolume*> &strong_volumes,
    Polygon &polygon,
    const PreparedPerimeter &prepared,
    const Layer *layer,
    const ModifierRegionsCache &slices_cache,
    PreciseSeamWarnings* warnings)
{
    assert(&prepared.polygon == &polygon);
    if (strong_volumes.empty() || layer == nullptr || !prepared.valid)
        return std::nullopt;
    const size_t layer_id = layer->id() - layer->object()->slicing_parameters().raft_layers();
    for (const ModelVolume *modifier : strong_volumes) {
        const auto it = slices_cache.find(modifier);
        if (it == slices_cache.end() || layer_id >= it->second.size())
            continue;
        const SegmentExtraction extracted = extract_perimeter_segments(
            prepared, it->second[layer_id], modifier->type(), {layer, modifier, warnings});
        // Full containment retains its existing skip policy, separately from segment selection.
        if (extracted.full_containment) {
            if (warnings)
                warnings->full_containment.store(true, std::memory_order_relaxed);
            continue;
        }
        // Policy: no usable segment passes the turn to the next modifier by priority. This includes a
        // modifier that crosses the perimeter but whose fragments were all discarded by binding: the
        // seam on that layer then comes from a lower-priority modifier (or weak/ordinary placement)
        // rather than from none, and the user is told through the "unable to process intersection"
        // warning. Deliberately not distinguished from "does not cross this perimeter".
        if (extracted.segments.empty())
            continue;
        if (warnings && extracted.segments.size() > 1)
            warnings->multiple_intersections.store(true, std::memory_order_relaxed);

        StrongSeamTarget target{Point(0, 0), 0};
        double longest = -1.;
        for (const PerimeterSegment &segment : extracted.segments) {
            assert(segment.strong_target.has_value());
            const StrongSeamTarget &candidate = *segment.strong_target;
            // Slice coordinates already include object rotation and retain bed axes.
            // Centering/instance translation cannot change rear (+Y), then left (-X) ordering.
            const bool farther_or_left = candidate.point.y() > target.point.y() ||
                (candidate.point.y() == target.point.y() && candidate.point.x() < target.point.x());
            if (segment.length > longest || (segment.length == longest && farther_or_left)) {
                longest = segment.length;
                target = candidate;
            }
        }
        auto result = insert_point_into_perimeter(target.point, target.edge_index, polygon);
        if (!result)
            return std::nullopt;
        // Preserve the insertion order: refining before first would shift the seam index.
        refine_at_vertex(result->second, +1, polygon);
        refine_at_vertex(result->second, -1, polygon);
        return result->first; // No later strong modifier or weak processing for this perimeter.
    }
    return std::nullopt;
}

// Convert ModelVolumeType of weak modifier to EnforcedBlockedSeamPoint.
// Precondition: called only with weak precise-seam types (filtered via is_precise_seam_weak()).
// Exhaustive switch (no default) so -Wswitch flags any future PRECISE_SEAM_* additions.
static EnforcedBlockedSeamPoint convert_weak_modifier_type(ModelVolumeType type) {
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

    std::vector<bool> segment_valid(result.size(), true);
    bool any_insertion_failed = false;

    for (const PointToInsert &pt : points_to_insert) {
        WeakModifierSegment &seg = result[pt.segment_idx];
        Point &point_coords = pt.is_left ? seg.left_point : seg.right_point;
        const size_t edge_idx = pt.position.edge_index;

        // Insert point with tolerance check
        std::optional<std::pair<Point, size_t>> insert_result =
            insert_point_into_perimeter(point_coords, edge_idx, polygon);

        if (insert_result.has_value()) {
            // Update coordinates in result (if point coincided with existing vertex, take its coordinates)
            point_coords = insert_result->first;
        } else {
            // Failed to insert point - segment becomes invalid
            segment_valid[pt.segment_idx] = false;
            any_insertion_failed = true;
        }
    }

    // Remove segments whose boundaries could not be inserted
    if (any_insertion_failed) {
        // Critical error: boundary point not inserted (shouldn't happen in normal conditions)
        BOOST_LOG_TRIVIAL(error) << "PreciseSeam: boundary point insertion failed, performing segment compaction";

        // Remove invalid segments (array compaction)
        size_t write_pos = 0;
        for (size_t read_pos = 0; read_pos < result.size(); ++read_pos) {
            if (segment_valid[read_pos]) {
                if (write_pos != read_pos) {
                    result[write_pos] = std::move(result[read_pos]);
                }
                ++write_pos;
            }
        }
        result.resize(write_pos);
    }

    // Group coincident boundaries by their snapped vertex, including wraparound to vertex 0.
    // Scan original vertices backwards so insertions cannot shift pending vertex indices.
    // O(vertices * segments), matching the boundary lookup below; typically only a few segments.
    for (size_t poly_idx = polygon.size(); poly_idx-- > 0; ) {
        bool refine_before = false;
        bool refine_after = false;
        for (const WeakModifierSegment &segment : result) {
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

    // Helper lambda: search for point index in modified polygon by coordinates.
    // Linear scan is intentional — O(N×M) is acceptable for typical M ≤ 5 weak segments.
    auto find_point_index = [&](const Point &pt) -> std::optional<size_t> {
        for (size_t i = 0; i < polygon.size(); ++i) {
            if (polygon[i] == pt) return i;
        }
        return std::nullopt;
    };

    // Apply types sequentially: segments are sorted low-priority-first
    // (bottom of object tree first), so higher-priority modifiers overwrite
    // lower-priority ones via last-write-wins.
    for (const auto &segment : result) {
        // Find boundary point indices in modified polygon
        std::optional<size_t> left_idx = find_point_index(segment.left_point);
        std::optional<size_t> right_idx = find_point_index(segment.right_point);

        if (!left_idx.has_value() || !right_idx.has_value()) {
            BOOST_LOG_TRIVIAL(error) << "PreciseSeam: boundary point not found in modified polygon, skipping segment";
            continue;
        }

        // Edge i joins vertices i and i+1, so the zone covers edges [left_idx, right_idx).
        // The edge starting at the right boundary lies outside the zone and must not be subdivided;
        // the 1 µm helper after the boundary keeps that edge too short to split, but this must not rely on it.
        // Edge types drive oversampling only: candidate types are assigned per point, including
        // both boundaries, in apply_weak_modifiers_to_perimeter. A zero-length zone marks no edge.
        for (size_t idx = left_idx.value(); idx != right_idx.value(); idx = (idx + 1) % polygon.size())
            edge_types[idx] = segment.type;
    }

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
    PreciseSeamWarnings* warnings)
{
    assert(&prepared.polygon == &polygon);
    std::vector<WeakModifierSegment> result;
    if (weak_volumes.empty() || layer == nullptr || !prepared.valid)
        return result;

    // Modifier cache indices exclude raft layers, unlike Layer::id().
    const size_t raft_layers = layer->object()->slicing_parameters().raft_layers();
    const size_t layer_id = layer->id() - raft_layers;
    for (const ModelVolume *modifier_volume : weak_volumes) {
        const auto it = slices_cache.find(modifier_volume);
        if (it == slices_cache.end() || layer_id >= it->second.size())
            continue;
        const ExtractionContext context{layer, modifier_volume, warnings};
        const SegmentExtraction extracted = extract_perimeter_segments(prepared, it->second[layer_id], modifier_volume->type(), context);
        // Keep the existing full-containment policy until it is changed explicitly.
        if (extracted.full_containment) {
            if (warnings)
                warnings->full_containment.store(true, std::memory_order_relaxed);
            continue;
        }
        const auto type = convert_weak_modifier_type(modifier_volume->type());
        for (const PerimeterSegment &segment : extracted.segments) {
            result.push_back({type, segment.polyline.points.front(), segment.begin,
                                   segment.polyline.points.back(), segment.end});
        }
    }
    // Boundary insertion, refinement and last-write-wins priority use the established path.
    return prepare_weak_modifier_segments(std::move(result), polygon);
}

// Apply weak modifier types to perimeter points based on segment boundaries.
// Find boundary points in refined polygon by coordinates and set types
// for all points inside each segment.
void apply_weak_modifiers_to_perimeter(
    const std::vector<WeakModifierSegment> &weak_segments,
    PrintObjectSeamData::LayerSeams &result,
    const SeamPlacerImpl::Perimeter &perimeter,
    bool &some_point_enforced)
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
