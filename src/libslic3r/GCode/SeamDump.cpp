// seamdump: the outer loops of a plate as SVG for an external seam choice. Design: docs/HLSD/precise-seam.md
#include "SeamDump.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <openssl/sha.h>
#include "libslic3r/AABBTreeLines.hpp"
#include "libslic3r/ArcFitter.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/IntersectionPoints.hpp"
#include "libslic3r/Line.hpp"
#include "libslic3r/MultiPoint.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"

namespace Slic3r {
namespace PreciseSeam {

namespace {

constexpr int    format_version = 1;
constexpr double tiny           = 1e-9; // mm; lengths below it are none.

int64_t to_hundredths(double mm) { return int64_t(std::llround(mm * 100.)); }
double  from_hundredths(int64_t h) { return double(h) / 100.; }

double cross2(const Vec2d &a, const Vec2d &b) { return a.x() * b.y() - a.y() * b.x(); }

std::string hex(const unsigned char *bytes, size_t size)
{
    static const char digits[] = "0123456789abcdef";
    std::string       out;
    for (size_t i = 0; i < size; ++i) {
        out.push_back(digits[bytes[i] >> 4]);
        out.push_back(digits[bytes[i] & 15]);
    }
    return out;
}

std::string sha256(const std::string &data)
{
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char *>(data.data()), data.size(), digest);
    return hex(digest, sizeof(digest));
}

std::string xml_escape(const std::string &text)
{
    std::string out;
    for (char c : text)
        switch (c) {
        case '&': out += "&amp;"; break;
        case '<': out += "&lt;"; break;
        case '>': out += "&gt;"; break;
        case '"': out += "&quot;"; break;
        case '\'': out += "&apos;"; break;
        default:
            // XML 1.0 has no other control characters than tab, CR and LF.
            if (static_cast<unsigned char>(c) >= 0x20 || c == '\t' || c == '\n' || c == '\r')
                out.push_back(c);
            else
                out.push_back(' ');
        }
    return out;
}

// ---------------------------------------------------------------------------------------------------------
// Canonical source loop: the published path is built from it, and its paint and length are measured on it.

struct CanonicalLoop {
    Points                    points; // Outer counterclockwise, hole clockwise; no repeated point.
    std::vector<Vec2d>        mm;     // The same points in mm.
    std::vector<double>       arc;    // Length from the first point to each point, the loop length last.
    std::vector<DumpPaintRun> paint;  // Within [0, length], split at the start.
    bool                      hole   = false;
    size_t                    source = 0;
    double                    area   = 0.; // Of the region inside, scaled units squared, positive.

    double length() const { return arc.back(); }
};

bool point_less(const Point &a, const Point &b) { return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y()); }

// Start of the lexicographically smallest rotation of a cyclic sequence, in linear time.
size_t least_rotation(const Points &points)
{
    const size_t n = points.size();
    size_t       i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        const Point &a = points[(i + k) % n], &b = points[(j + k) % n];
        if (a == b) {
            ++k;
            continue;
        }
        if (point_less(b, a))
            i += k + 1;
        else
            j += k + 1;
        if (i == j)
            ++j;
        k = 0;
    }
    return std::min(i, j);
}

// Paint runs of a loop of `length` moved to positions s -> s' and kept within [0, length], split at the start.
std::vector<DumpPaintRun> remap_paint(const std::vector<DumpPaintRun> &runs, double length, bool reverse, double shift)
{
    std::vector<DumpPaintRun> out;
    if (length <= tiny)
        return out;
    for (DumpPaintRun run : runs) {
        if (run.end - run.begin >= length - tiny) {
            out.push_back({0., length, run.blocked});
            continue;
        }
        if (reverse)
            run = {length - run.end, length - run.begin, run.blocked};
        double begin = std::fmod(run.begin - shift, length);
        if (begin < 0.)
            begin += length;
        const double end = begin + (run.end - run.begin);
        if (end <= length)
            out.push_back({begin, end, run.blocked});
        else {
            out.push_back({begin, length, run.blocked});
            out.push_back({0., end - length, run.blocked});
        }
    }
    return out;
}

std::optional<CanonicalLoop> canonicalize(const DumpLoop &loop, size_t source)
{
    CanonicalLoop out;
    out.hole   = loop.hole;
    out.source = source;
    // Repeated points: the shared ends of paths and a closing point. They add no length, so paint stays put.
    for (const Point &p : loop.polygon.points)
        if (out.points.empty() || out.points.back() != p)
            out.points.push_back(p);
    while (out.points.size() > 1 && out.points.back() == out.points.front())
        out.points.pop_back();
    if (out.points.size() < 3)
        return std::nullopt;
    const double area = Polygon(out.points).area();
    if (area == 0.)
        return std::nullopt;
    std::vector<double> source_arc{0.};
    for (size_t i = 0; i < out.points.size(); ++i)
        source_arc.push_back(source_arc.back() + (unscale(out.points[(i + 1) % out.points.size()]) - unscale(out.points[i])).norm());
    const double length = source_arc.back();
    // Outer loops counterclockwise, holes clockwise, whatever the print direction.
    const bool reverse = (area > 0.) == loop.hole;
    if (reverse)
        std::reverse(out.points.begin() + 1, out.points.end());
    const size_t start = least_rotation(out.points);
    // Length from the first point to the new start, along the (possibly reversed) loop.
    double shift = 0.;
    for (size_t i = 0; i < start; ++i)
        shift += (unscale(out.points[i + 1]) - unscale(out.points[i])).norm();
    std::rotate(out.points.begin(), out.points.begin() + start, out.points.end());
    for (const Point &p : out.points)
        out.mm.push_back(unscale(p));
    out.arc = {0.};
    for (size_t i = 0; i < out.mm.size(); ++i)
        out.arc.push_back(out.arc.back() + (out.mm[(i + 1) % out.mm.size()] - out.mm[i]).norm());
    out.paint = remap_paint(loop.paint, length, reverse, shift);
    out.area  = std::abs(area);
    return out;
}

// ---------------------------------------------------------------------------------------------------------
// Published path: what the file writes, on the 0.01 mm grid. `d` is its only definition.

struct Segment {
    Vec2i64 a, b;          // Hundredths of a mm.
    bool    arc   = false;
    int64_t r     = 0;     // Arc radius in hundredths.
    bool    ccw   = false; // Arc direction, Y up.
    bool    large = false; // Arc angle above pi.
    Vec2d   center{0., 0.}; // mm, derived from a, b, r and the flags as SVG does.
    double  angle = 0.;    // Arc angle in radians, from the derived center.

    Vec2d start() const { return Vec2d(from_hundredths(a.x()), from_hundredths(a.y())); }
    Vec2d end() const { return Vec2d(from_hundredths(b.x()), from_hundredths(b.y())); }
    double length() const { return arc ? from_hundredths(r) * angle : (end() - start()).norm(); }
    // Unit direction at the start or the end.
    Vec2d tangent(bool at_end) const
    {
        if (!arc)
            return (end() - start()).normalized();
        const Vec2d radial = (at_end ? end() : start()) - center;
        const Vec2d t      = ccw ? Vec2d(-radial.y(), radial.x()) : Vec2d(radial.y(), -radial.x());
        return t.normalized();
    }
    Vec2d point_at(double s) const
    {
        const double l = length();
        if (!arc)
            return l > tiny ? Vec2d(start() + (end() - start()) * (s / l)) : start();
        const Vec2d  radial = start() - center;
        const double turn   = (ccw ? 1. : -1.) * (l > tiny ? angle * s / l : 0.);
        return center + Vec2d(radial.x() * std::cos(turn) - radial.y() * std::sin(turn), radial.x() * std::sin(turn) + radial.y() * std::cos(turn));
    }
};

// Sets the center and angle of an arc from its rounded ends, radius and flags, as an SVG reader does. False if
// SVG would have to enlarge the radius: then the arc is no arc.
bool derive_arc(Segment &s)
{
    const Vec2d  a = s.start(), b = s.end();
    const double r = from_hundredths(s.r);
    const double half = (b - a).norm() / 2.;
    if (half <= tiny || r < half)
        return false;
    const Vec2d  mid  = (a + b) / 2.;
    const Vec2d  dir  = (b - a).normalized();
    const Vec2d  left = Vec2d(-dir.y(), dir.x());
    const double h    = std::sqrt(std::max(0., r * r - half * half));
    // A small counterclockwise arc has its center left of the chord, a large one right of it; clockwise mirrors.
    s.center = mid + left * ((s.ccw != s.large) ? h : -h);
    const double from = std::atan2(a.y() - s.center.y(), a.x() - s.center.x());
    const double to   = std::atan2(b.y() - s.center.y(), b.x() - s.center.x());
    double       turn = s.ccw ? to - from : from - to;
    while (turn <= 0.)
        turn += 2. * PI;
    while (turn > 2. * PI)
        turn -= 2. * PI;
    s.angle = turn;
    return true;
}

struct Published {
    std::vector<Segment> segments;
    std::vector<double>  arc; // Length from M to the start of each segment, the path length last.

    double length() const { return arc.back(); }
    Vec2d  point_at(double s) const
    {
        s = std::clamp(s, 0., length());
        const size_t i = std::min(size_t(std::upper_bound(arc.begin(), arc.end(), s) - arc.begin()), segments.size()) - 1;
        return segments[i].point_at(s - arc[i]);
    }
    // Arc length of the point of the path nearest to `p`.
    double project(const Vec2d &p) const
    {
        double best = std::numeric_limits<double>::max(), at = 0.;
        for (size_t i = 0; i < segments.size(); ++i) {
            const Segment &seg = segments[i];
            double         s   = 0.;
            if (!seg.arc) {
                const Vec2d  ab = seg.end() - seg.start();
                const double l2 = ab.squaredNorm();
                s               = l2 > 0. ? std::clamp((p - seg.start()).dot(ab) / l2, 0., 1.) * std::sqrt(l2) : 0.;
            } else {
                // The angle of p from the start, in the arc direction, clamped to the arc.
                const Vec2d  u = seg.start() - seg.center, v = p - seg.center;
                double       t = std::atan2(cross2(u, v), u.dot(v));
                if (!seg.ccw)
                    t = -t;
                if (t < 0.)
                    t += 2. * PI;
                // Beyond the arc: the nearer end.
                if (t > seg.angle)
                    t = (t - seg.angle < 2. * PI - t) ? seg.angle : 0.;
                s = from_hundredths(seg.r) * t;
            }
            const double d = (seg.point_at(s) - p).squaredNorm();
            if (d < best) {
                best = d;
                at   = arc[i] + s;
            }
        }
        return at;
    }
    // The path as a polygon in mm: arcs split so that no chord leaves them by more than 0.5 um.
    std::vector<Vec2d> polygonize() const
    {
        std::vector<Vec2d> out;
        for (const Segment &seg : segments) {
            out.push_back(seg.start());
            if (seg.arc) {
                const double r    = from_hundredths(seg.r);
                const double step = 2. * std::acos(std::max(0., 1. - 0.0005 / r));
                const size_t n    = std::max<size_t>(2, size_t(std::ceil(seg.angle / std::max(step, 1e-3))));
                for (size_t k = 1; k < n; ++k)
                    out.push_back(seg.point_at(seg.length() * double(k) / double(n)));
            }
        }
        return out;
    }
};

enum class Simplification { Arcs, Lines, None };

// The published path of a loop: simplified with T/2, then put on the grid. Empty if nothing is left.
std::optional<Published> publish(const CanonicalLoop &loop, Simplification mode, double tolerance)
{
    Points closed = loop.points;
    closed.push_back(loop.points.front());
    std::vector<PathFittingData> fitting;
    if (mode == Simplification::Arcs)
        ArcFitter::do_arc_fitting_and_simplify(closed, fitting, scale_(tolerance / 2.));
    else {
        if (mode == Simplification::Lines)
            closed = MultiPoint::_douglas_peucker(closed, scale_(tolerance / 2.));
        fitting.push_back({0, closed.size() - 1, EMovePathType::Linear_move, ArcSegment()});
    }
    const auto grid = [](const Point &p) { return Vec2i64(to_hundredths(unscale<double>(p.x())), to_hundredths(unscale<double>(p.y()))); };
    std::vector<Segment> segments;
    const auto add_line = [&segments](const Vec2i64 &a, const Vec2i64 &b) {
        // Ends that meet on the grid are one point.
        if (a != b)
            segments.push_back({a, b});
    };
    for (PathFittingData &part : fitting) {
        if (part.is_arc_move() && part.arc_data.is_valid()) {
            Segment s;
            s.a     = grid(closed[part.start_point_index]);
            s.b     = grid(closed[part.end_point_index]);
            s.arc   = true;
            s.r     = to_hundredths(unscale<double>(part.arc_data.radius));
            s.ccw   = part.arc_data.direction == ArcDirection::Arc_Dir_CCW;
            s.large = std::abs(part.arc_data.angle_radians) > PI;
            // An arc whose ends meet on the grid, or whose radius SVG would enlarge, becomes a line.
            if (s.a != s.b && derive_arc(s))
                segments.push_back(s);
            else
                add_line(s.a, s.b);
        } else
            for (size_t i = part.start_point_index; i < part.end_point_index; ++i)
                add_line(grid(closed[i]), grid(closed[i + 1]));
    }
    if (segments.size() < 2)
        return std::nullopt;
    Published out;
    out.segments = std::move(segments);
    out.arc      = {0.};
    for (const Segment &s : out.segments)
        out.arc.push_back(out.arc.back() + s.length());
    return out;
}

// Signed area of a polygon in mm, positive counterclockwise.
double signed_area(const std::vector<Vec2d> &polygon)
{
    double area = 0.;
    for (size_t i = 0; i < polygon.size(); ++i)
        area += cross2(polygon[i], polygon[(i + 1) % polygon.size()]);
    return area / 2.;
}

Polygon scaled_polygon(const std::vector<Vec2d> &polygon)
{
    Points points;
    for (const Vec2d &p : polygon)
        points.push_back(Point::new_scale(p));
    return Polygon(points);
}

std::vector<Linef> closed_lines(const std::vector<Vec2d> &polygon)
{
    std::vector<Linef> lines;
    for (size_t i = 0; i < polygon.size(); ++i)
        lines.emplace_back(polygon[i], polygon[(i + 1) % polygon.size()]);
    return lines;
}

// Points along a closed polygon, at most `step` apart, its vertices included.
std::vector<Vec2d> sample(const std::vector<Vec2d> &polygon, double step)
{
    std::vector<Vec2d> out;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const Vec2d &a = polygon[i], &b = polygon[(i + 1) % polygon.size()];
        const size_t n = std::max<size_t>(1, size_t(std::ceil((b - a).norm() / step)));
        for (size_t k = 0; k < n; ++k)
            out.push_back(a + (b - a) * (double(k) / double(n)));
    }
    return out;
}

// Whether the published path keeps within T of the source loop both ways, sampled at T/4 on either side, and
// keeps the orientation of the loop with an area.
bool within_tolerance(const CanonicalLoop &loop, const std::vector<Vec2d> &published, double tolerance)
{
    const double area = signed_area(published);
    if (std::abs(area) <= tiny || (area > 0.) == loop.hole)
        return false;
    const AABBTreeLines::LinesDistancer<Linef> to_source(closed_lines(loop.mm)), to_published(closed_lines(published));
    const double step = tolerance / 4.;
    for (const Vec2d &p : sample(published, step))
        if (to_source.distance_from_lines<false>(p) > tolerance)
            return false;
    for (const Vec2d &p : sample(loop.mm, step))
        if (to_published.distance_from_lines<false>(p) > tolerance)
            return false;
    return true;
}

// ---------------------------------------------------------------------------------------------------------
// Layers: loops of one layer with their published form.

struct LayerLoop {
    CanonicalLoop            loop;
    Simplification           mode = Simplification::Arcs;
    std::optional<Published> path;     // None once dropped.
    std::vector<Vec2d>       polygon;  // Polygonized path.
    int64_t                  id = 0;   // L number.
    std::vector<int64_t>     from;     // Former L numbers after a split or merge.
    std::optional<size_t>    parent;   // Index of the outer loop around a hole, in this layer.
};

struct BuiltLayer {
    const DumpLayer          *source = nullptr;
    std::vector<LayerLoop>    loops;
    std::vector<std::string>  dropped; // "S<n>:reason" before numbering is known.
};

// Builds the published path of each loop with the mildest simplification that keeps the tolerance and the
// topology of the layer: no crossing between published loops and each hole inside its outer loop.
void publish_layer(BuiltLayer &layer, double tolerance)
{
    const auto try_mode = [&](LayerLoop &l) {
        for (;;) {
            l.path = publish(l.loop, l.mode, tolerance);
            if (l.path) {
                l.polygon = l.path->polygonize();
                if (within_tolerance(l.loop, l.polygon, tolerance))
                    return;
            }
            if (l.mode == Simplification::None) {
                l.path.reset();
                return;
            }
            l.mode = Simplification(int(l.mode) + 1);
        }
    };
    for (LayerLoop &l : layer.loops) {
        try_mode(l);
        if (!l.path)
            layer.dropped.push_back("S" + std::to_string(l.loop.source) + ":tolerance");
    }
    // Crossings and nesting: each loop involved falls back to milder simplification, until none is left.
    for (;;) {
        std::vector<bool> bad(layer.loops.size(), false);
        Lines             lines;
        std::vector<size_t> owner;
        // A far, empty line between loops: the intersection search skips neighbouring lines.
        const Point separator(std::numeric_limits<coord_t>::max() / 4, std::numeric_limits<coord_t>::max() / 4);
        for (size_t i = 0; i < layer.loops.size(); ++i) {
            const LayerLoop &l = layer.loops[i];
            if (!l.path)
                continue;
            for (size_t k = 0; k < l.polygon.size(); ++k) {
                lines.emplace_back(Point::new_scale(l.polygon[k]), Point::new_scale(l.polygon[(k + 1) % l.polygon.size()]));
                owner.push_back(i);
            }
            lines.emplace_back(separator, separator);
            owner.push_back(size_t(-1));
        }
        for (const IntersectionLines &x : get_intersections(lines))
            for (uint32_t li : {x.line_index1, x.line_index2})
                if (owner[li] != size_t(-1))
                    bad[owner[li]] = true;
        for (size_t i = 0; i < layer.loops.size(); ++i) {
            const LayerLoop &l = layer.loops[i];
            if (l.path && l.parent && layer.loops[*l.parent].path &&
                !scaled_polygon(layer.loops[*l.parent].polygon).contains(Point::new_scale(l.polygon.front())))
                bad[i] = true;
        }
        bool changed = false;
        for (size_t i = 0; i < layer.loops.size(); ++i)
            if (bad[i]) {
                LayerLoop &l = layer.loops[i];
                if (l.mode == Simplification::None) {
                    l.path.reset();
                    layer.dropped.push_back("S" + std::to_string(l.loop.source) + ":topology");
                } else {
                    l.mode = Simplification(int(l.mode) + 1);
                    try_mode(l);
                    if (!l.path)
                        layer.dropped.push_back("S" + std::to_string(l.loop.source) + ":tolerance");
                }
                changed = true;
            }
        if (!changed)
            break;
    }
}

// ---------------------------------------------------------------------------------------------------------
// Numbers across layers.

// Area of the region inside a canonical loop, as a counterclockwise polygon.
Polygon region(const CanonicalLoop &loop)
{
    Polygon p(loop.points);
    if (loop.hole)
        p.reverse();
    return p;
}

bool overlaps(const CanonicalLoop &a, const CanonicalLoop &b)
{
    const Polygon pa = region(a), pb = region(b);
    if (!get_extents(pa).overlap(get_extents(pb)))
        return false;
    double area = 0.;
    for (const Polygon &p : intersection(Polygons{pa}, Polygons{pb}))
        area += p.area();
    return area > 0.;
}

bool sequence_less(const Points &a, const Points &b)
{
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end(), point_less);
}

// L numbers of one object: a loop keeps the number of the only loop of its kind it overlaps on the layer below,
// when that loop overlaps no other; every other loop gets a new number, with the loops it came from.
void number_loops(std::vector<BuiltLayer> &layers)
{
    int64_t next = 1;
    for (size_t li = 0; li < layers.size(); ++li) {
        std::vector<LayerLoop> &current = layers[li].loops;
        std::vector<LayerLoop> *below   = li > 0 ? &layers[li - 1].loops : nullptr;
        // Components of the overlap graph per kind, as a union-find over below + current.
        const size_t       nb = below ? below->size() : 0;
        std::vector<size_t> parent(nb + current.size());
        std::iota(parent.begin(), parent.end(), 0);
        const auto find = [&parent](size_t i) {
            while (parent[i] != i)
                i = parent[i] = parent[parent[i]];
            return i;
        };
        std::vector<bool> linked(nb + current.size(), false);
        for (size_t c = 0; c < current.size(); ++c)
            for (size_t b = 0; b < nb; ++b)
                if ((*below)[b].loop.hole == current[c].loop.hole && overlaps((*below)[b].loop, current[c].loop)) {
                    parent[find(nb + c)] = find(b);
                    linked[nb + c] = linked[b] = true;
                }
        std::map<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>> components; // root -> (below, current)
        for (size_t b = 0; b < nb; ++b)
            if (linked[b])
                components[find(b)].first.push_back(b);
        for (size_t c = 0; c < current.size(); ++c)
            if (linked[nb + c])
                components[find(nb + c)].second.push_back(c);
        std::vector<size_t> fresh;
        for (size_t c = 0; c < current.size(); ++c)
            if (!linked[nb + c])
                fresh.push_back(c);
        for (auto &[root, members] : components) {
            if (members.first.size() == 1 && members.second.size() == 1) {
                current[members.second.front()].id = (*below)[members.first.front()].id;
                continue;
            }
            std::vector<int64_t> from;
            for (size_t b : members.first)
                from.push_back((*below)[b].id);
            std::sort(from.begin(), from.end());
            for (size_t c : members.second) {
                current[c].from = from;
                fresh.push_back(c);
            }
        }
        // New numbers in a fixed order: outer loops first, then by their points, then by source.
        std::sort(fresh.begin(), fresh.end(), [&current](size_t a, size_t b) {
            const CanonicalLoop &la = current[a].loop, &lb = current[b].loop;
            if (la.hole != lb.hole)
                return !la.hole;
            if (la.points != lb.points)
                return sequence_less(la.points, lb.points);
            return la.source < lb.source;
        });
        for (size_t c : fresh)
            current[c].id = next++;
    }
}

// The outer loop around each hole: the smallest one of the layer that contains its first point.
void find_parents(BuiltLayer &layer)
{
    for (LayerLoop &l : layer.loops) {
        if (!l.loop.hole)
            continue;
        double best = std::numeric_limits<double>::max();
        for (size_t i = 0; i < layer.loops.size(); ++i) {
            const CanonicalLoop &outer = layer.loops[i].loop;
            if (!outer.hole && outer.area < best && Polygon(outer.points).contains(l.loop.points.front())) {
                best     = outer.area;
                l.parent = i;
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------------------
// Candidates and paint on the published path.

struct Candidate {
    double s = 0.;
    Vec2d  at;
    int    angle  = 180; // Degrees on the side of the material.
    bool   convex = false;
    double turn   = 0.;  // Signed turn, left positive.
};

// Every vertex of the path with a turn, ranked by its turn and then its distance from M; the first `count` of
// them in path order. Material lies left of the path for outer loops and holes alike.
std::vector<Candidate> candidates(const Published &path, int count)
{
    std::vector<Candidate> all;
    const size_t           n = path.segments.size();
    for (size_t i = 0; i < n; ++i) {
        const Vec2d  in   = path.segments[(i + n - 1) % n].tangent(true);
        const Vec2d  out  = path.segments[i].tangent(false);
        const double turn = std::atan2(cross2(in, out), in.dot(out));
        Candidate    c;
        c.s      = path.arc[i];
        c.at     = path.segments[i].start();
        c.turn   = turn;
        c.angle  = int(std::lround(180. - turn * 180. / PI));
        c.convex = turn > 0.;
        // A vertex where the path runs straight on is no corner.
        if (c.angle != 180)
            all.push_back(c);
    }
    std::stable_sort(all.begin(), all.end(), [](const Candidate &a, const Candidate &b) {
        return std::abs(a.turn) != std::abs(b.turn) ? std::abs(a.turn) > std::abs(b.turn) : a.s < b.s;
    });
    if (all.size() > size_t(std::max(count, 0)))
        all.resize(size_t(std::max(count, 0)));
    std::sort(all.begin(), all.end(), [](const Candidate &a, const Candidate &b) { return a.s < b.s; });
    return all;
}

struct PaintInterval {
    int64_t a = 0, b = 0; // Hundredths along the published path, a < b.
    bool    blocked = false;
};

// 0 none, 1 enforced, 2 blocked at a position of the canonical loop.
int paint_state(const std::vector<DumpPaintRun> &runs, double s)
{
    int state = 0;
    for (const DumpPaintRun &run : runs)
        if (s >= run.begin && s <= run.end)
            state = std::max(state, run.blocked ? 2 : 1);
    return state;
}

// Arc length on the canonical loop of the point nearest to `p`.
double project_source(const CanonicalLoop &loop, const Vec2d &p)
{
    double best = std::numeric_limits<double>::max(), at = 0.;
    for (size_t i = 0; i < loop.mm.size(); ++i) {
        const Vec2d  a = loop.mm[i], ab = loop.mm[(i + 1) % loop.mm.size()] - a;
        const double l2 = ab.squaredNorm();
        const double t  = l2 > 0. ? std::clamp((p - a).dot(ab) / l2, 0., 1.) : 0.;
        const double d  = (a + ab * t - p).squaredNorm();
        if (d < best) {
            best = d;
            at   = loop.arc[i] + t * (loop.arc[i + 1] - loop.arc[i]);
        }
    }
    return at;
}

// Paint moved to the published path: each end of a run goes to the nearest point of the path, and each piece
// between two ends takes the paint of the source point nearest to its middle. An approximation: 0.01 mm is the
// step of the written ends, not a bound of their error.
std::vector<PaintInterval> publish_paint(const CanonicalLoop &loop, const Published &path)
{
    std::vector<PaintInterval> out;
    if (loop.paint.empty())
        return out;
    const double        length = path.length();
    std::vector<double> cuts{0., length};
    for (const DumpPaintRun &run : loop.paint)
        for (double s : {run.begin, run.end}) {
            const double at = std::fmod(s, loop.length());
            const size_t i  = std::min(size_t(std::upper_bound(loop.arc.begin(), loop.arc.end(), at) - loop.arc.begin()), loop.mm.size()) - 1;
            const double l  = loop.arc[i + 1] - loop.arc[i];
            const Vec2d  p  = loop.mm[i] + (loop.mm[(i + 1) % loop.mm.size()] - loop.mm[i]) * (l > 0. ? (at - loop.arc[i]) / l : 0.);
            cuts.push_back(path.project(p));
        }
    std::sort(cuts.begin(), cuts.end());
    for (size_t i = 0; i + 1 < cuts.size(); ++i) {
        if (cuts[i + 1] - cuts[i] <= tiny)
            continue;
        const int state = paint_state(loop.paint, project_source(loop, path.point_at((cuts[i] + cuts[i + 1]) / 2.)));
        if (state == 0)
            continue;
        const int64_t a = to_hundredths(cuts[i]), b = to_hundredths(cuts[i + 1]);
        if (!out.empty() && out.back().b == a && out.back().blocked == (state == 2))
            out.back().b = b;
        else if (b > a)
            out.push_back({a, b, state == 2});
    }
    return out;
}

// ---------------------------------------------------------------------------------------------------------
// Writing.

std::string number(double mm) { return seamdump_number(to_hundredths(mm)); }
std::string point_text(const Vec2i64 &p) { return seamdump_number(p.x()) + " " + seamdump_number(p.y()); }

struct LoopText {
    std::string path;  // The <path> element of the loop.
    std::string paint; // Helper paths that draw its paint for a person.
};

LoopText loop_text(const std::string &object_id, size_t layer_number, const LayerLoop &l, const BuiltLayer &layer, int count)
{
    const Published &path = *l.path;
    std::string      d    = "M" + point_text(path.segments.front().a);
    std::string      arcs;
    for (size_t i = 0; i < path.segments.size(); ++i) {
        const Segment &s = path.segments[i];
        // The last segment ends at M: a line closes with Z alone, an arc is written and then closed.
        const bool last = i + 1 == path.segments.size();
        if (s.arc) {
            d += " A" + seamdump_number(s.r) + " " + seamdump_number(s.r) + " 0 " + (s.large ? "1" : "0") + " " + (s.ccw ? "1" : "0") + " " + point_text(s.b);
            if (!arcs.empty())
                arcs += ";";
            arcs += number(s.center.x()) + "," + number(s.center.y()) + "," + seamdump_number(s.r);
        } else if (!last)
            d += " L" + point_text(s.b);
    }
    d += " Z";
    const std::vector<PaintInterval> paint = publish_paint(l.loop, path);
    std::string                      c;
    for (const Candidate &cand : candidates(path, count)) {
        if (!c.empty())
            c += ";";
        const int64_t s    = to_hundredths(cand.s);
        char          mark = 0;
        for (const PaintInterval &p : paint)
            if (s >= p.a && s <= p.b && (mark != 'b'))
                mark = p.blocked ? 'b' : 'e';
        c += number(cand.at.x()) + "," + number(cand.at.y()) + "," + std::to_string(cand.angle) + "," + (cand.convex ? "1" : "0");
        if (mark)
            c += std::string(",") + mark;
    }
    std::string paint_attr;
    for (const PaintInterval &p : paint)
        paint_attr += std::string(paint_attr.empty() ? "" : ";") + (p.blocked ? "block:" : "enforce:") + seamdump_number(p.a) + "-" + seamdump_number(p.b);

    const std::string id = object_id + "-Z" + std::to_string(layer_number) + "-L" + std::to_string(l.id);
    LoopText          out;
    out.path = "<path id=\"" + id + "\" data-source-loop=\"" + std::to_string(l.loop.source) + "\" data-kind=\"" + (l.loop.hole ? "hole" : "outer") + "\"";
    if (l.parent && layer.loops[*l.parent].path)
        out.path += " data-parent=\"L" + std::to_string(layer.loops[*l.parent].id) + "\"";
    if (!l.from.empty()) {
        out.path += " data-from=\"";
        for (size_t i = 0; i < l.from.size(); ++i)
            out.path += (i ? " L" : "L") + std::to_string(l.from[i]);
        out.path += "\"";
    }
    out.path += " data-len=\"" + number(path.length()) + "\"";
    if (!arcs.empty())
        out.path += " data-arcs=\"" + arcs + "\"";
    if (!c.empty())
        out.path += " data-c=\"" + c + "\"";
    if (!paint_attr.empty())
        out.path += " data-paint=\"" + paint_attr + "\"";
    out.path += " d=\"" + d + "\"/>\n";
    // For a person: the painted pieces in colour, as lines through the path; no part of the loop data.
    for (const PaintInterval &p : paint) {
        std::string piece;
        const double a = from_hundredths(p.a), b = from_hundredths(p.b);
        const size_t n = std::max<size_t>(1, size_t(std::ceil((b - a) / 0.5)));
        for (size_t k = 0; k <= n; ++k) {
            const Vec2d q = path.point_at(a + (b - a) * double(k) / double(n));
            piece += (k ? " L" : "M") + number(q.x()) + " " + number(q.y());
        }
        out.paint += "<path class=\"" + std::string(p.blocked ? "block" : "enforce") + "\" data-loop=\"L" + std::to_string(l.id) +
                     "\" stroke=\"" + (p.blocked ? "red" : "green") + "\" stroke-width=\"0.3\" d=\"" + piece + "\"/>\n";
    }
    return out;
}

// Source records of the objects, unrounded, for the fingerprint: a later import compares them with its slicing.
std::string fingerprint_text(const std::vector<DumpObject> &objects, const DumpSettings &settings)
{
    const auto nm = [](double mm) { return std::to_string(std::llround(mm * 1e6)); };
    std::string text = "orca-seamdump-fingerprint " + std::to_string(format_version) + "\nreach " + nm(settings.reach) + "\n";
    for (size_t o = 0; o < objects.size(); ++o) {
        const DumpObject &object = objects[o];
        text += "object " + std::to_string(o + 1) + " " + object.name + "\n";
        for (const Vec2d &shift : object.instances)
            text += "instance " + nm(shift.x()) + " " + nm(shift.y()) + "\n";
        for (const DumpLayer &layer : object.layers) {
            text += "layer " + std::to_string(layer.number) + " " + nm(layer.slice_z) + " " + nm(layer.print_z) + " " + nm(layer.height) + "\n";
            for (size_t l = 0; l < layer.loops.size(); ++l) {
                const DumpLoop &loop = layer.loops[l];
                text += "loop " + std::to_string(l) + (loop.hole ? " hole" : " outer");
                for (const Point &p : loop.polygon.points)
                    text += " " + std::to_string(p.x()) + "," + std::to_string(p.y());
                text += "\n";
                for (const DumpPaintRun &run : loop.paint)
                    text += std::string("paint ") + (run.blocked ? "block " : "enforce ") + nm(run.begin) + " " + nm(run.end) + "\n";
            }
        }
    }
    return text;
}

struct Built {
    std::string body;    // Object groups.
    BoundingBoxf bounds; // Of the published loops, mm.
    size_t      dropped = 0;
    std::string dropped_list;
};

Built build(const std::vector<DumpObject> &objects, double tolerance, int count)
{
    Built out;
    for (size_t o = 0; o < objects.size(); ++o) {
        const DumpObject &object    = objects[o];
        const std::string object_id = "O" + std::to_string(o + 1);
        std::vector<BuiltLayer> layers(object.layers.size());
        for (size_t li = 0; li < object.layers.size(); ++li) {
            BuiltLayer &layer = layers[li];
            layer.source      = &object.layers[li];
            for (size_t s = 0; s < layer.source->loops.size(); ++s) {
                if (std::optional<CanonicalLoop> canonical = canonicalize(layer.source->loops[s], s))
                    layer.loops.push_back({std::move(*canonical)});
                else
                    layer.dropped.push_back("S" + std::to_string(s) + ":degenerate");
            }
            find_parents(layer);
            publish_layer(layer, tolerance);
        }
        number_loops(layers);

        std::string instances;
        for (const Vec2d &shift : object.instances)
            instances += (instances.empty() ? "" : ";") + number(shift.x()) + "," + number(shift.y());
        out.body += " <g id=\"" + object_id + "\" data-object=\"" + xml_escape(object.name) + "\" data-instances=\"" + instances + "\">\n";
        std::string last_record, last_id; // The last layer written with geometry.
        for (const BuiltLayer &layer : layers) {
            const size_t number_ = layer.source->number;
            for (const std::string &d : layer.dropped) {
                out.dropped_list += (out.dropped_list.empty() ? "" : " ") + object_id + "-Z" + std::to_string(number_) + "-" + d;
                ++out.dropped;
            }
            std::string loops, paint, record;
            for (const LayerLoop &l : layer.loops) {
                if (!l.path)
                    continue;
                const LoopText text = loop_text(object_id, number_, l, layer, count);
                loops += "   " + text.path;
                paint += text.paint;
                for (const Vec2d &p : l.polygon)
                    out.bounds.merge(p);
            }
            if (loops.empty())
                continue;
            // The record of a layer without its own numbers: equal records repeat the same loops.
            const std::string layer_id = object_id + "-Z" + std::to_string(number_);
            record                     = loops;
            const std::string prefix   = layer_id + "-L";
            for (size_t pos = 0; (pos = record.find(prefix, pos)) != std::string::npos;)
                record.replace(pos, layer_id.size(), "#");
            const std::string head = "  <g id=\"" + layer_id + "\" data-layer=\"" + std::to_string(number_) + "\" data-z=\"" + number(layer.source->print_z) + "\"";
            if (!last_id.empty() && record == last_record) {
                out.body += head + " data-same=\"" + last_id + "\"/>\n";
                continue;
            }
            out.body += head + ">\n" + loops;
            if (!paint.empty())
                out.body += "   " + paint;
            out.body += "  </g>\n";
            last_record = record;
            last_id     = layer_id;
        }
        out.body += " </g>\n";
    }
    return out;
}

std::string replace_all(std::string text, const std::string &from, const std::string &to)
{
    for (size_t pos = 0; (pos = text.find(from, pos)) != std::string::npos; pos += to.size())
        text.replace(pos, from.size(), to);
    return text;
}

// The whole file: its size and token estimate count themselves, then the dump id hashes the bytes with zeros
// in its place. The placeholders hold control characters, which xml_escape() never leaves in a name.
std::string assemble(const Built &built, const std::string &fingerprint, double tolerance, int count)
{
    BoundingBoxf box = built.bounds;
    if (!box.defined)
        box = BoundingBoxf(Vec2d(0., 0.), Vec2d(0., 0.));
    const double x0 = std::floor(box.min.x()) - 1., y0 = std::floor(box.min.y()) - 1.;
    const double x1 = std::ceil(box.max.x()) + 1., y1 = std::ceil(box.max.y()) + 1.;
    const std::string zeros(64, '0');
    std::string       head = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"" + number(x0) + " " + number(y0) + " " + number(x1 - x0) + " " +
                       number(y1 - y0) + "\" width=\"" + number(x1 - x0) + "mm\" height=\"" + number(y1 - y0) + "mm\"" +
                       " data-format=\"orca-seamdump\" data-version=\"" + std::to_string(format_version) + "\" data-dump-id=\"" + zeros +
                       "\" data-fingerprint=\"" + fingerprint + "\" data-tolerance=\"" + number(tolerance) + "\" data-candidates-max=\"" +
                       std::to_string(count) + "\" data-layer-origin=\"model\" data-bytes=\"\x01" "B" "\x01\" data-token-estimate=\"\x01" "E" "\x01\" data-dropped=\"" +
                       std::to_string(built.dropped) + "\"" + (built.dropped_list.empty() ? "" : " data-dropped-loops=\"" + built.dropped_list + "\"") + ">\n";
    // Y up, as on the bed: the group is mirrored within the same view box.
    const std::string body = head + " <g transform=\"matrix(1 0 0 -1 0 " + number(y0 + y1) + ")\" fill=\"none\" stroke=\"black\" stroke-width=\"0.1\">\n" +
                             built.body + " </g>\n</svg>\n";
    // The size field is part of the size: settle it, then the estimate, which cannot change the size by more.
    std::string text;
    size_t      bytes = body.size();
    for (int i = 0; i < 8; ++i) {
        text = replace_all(replace_all(body, "\x01" "B" "\x01", std::to_string(bytes)), "\x01" "E" "\x01", std::to_string(bytes / 4));
        if (text.size() == bytes)
            break;
        bytes = text.size();
    }
    const std::string id = sha256(text);
    return replace_all(text, "data-dump-id=\"" + zeros + "\"", "data-dump-id=\"" + id + "\"");
}

} // namespace

std::string seamdump_number(int64_t hundredths)
{
    const bool     negative = hundredths < 0;
    const uint64_t a        = negative ? uint64_t(-(hundredths + 1)) + 1 : uint64_t(hundredths);
    std::string    out      = (negative ? "-" : "") + std::to_string(a / 100);
    const int      fraction = int(a % 100);
    if (fraction != 0) {
        out += '.';
        out += char('0' + fraction / 10);
        if (fraction % 10 != 0)
            out += char('0' + fraction % 10);
    }
    return out;
}

DumpResult write_seamdump(const std::vector<DumpObject> &objects, const DumpSettings &settings)
{
    const std::string fingerprint = sha256(fingerprint_text(objects, settings));
    // Above the limit: 2T, 4T, then half of K once; K never grows.
    std::vector<std::pair<double, int>> steps{{settings.tolerance, settings.candidates}};
    if (settings.limit != 0) {
        steps.push_back({settings.tolerance * 2., settings.candidates});
        steps.push_back({settings.tolerance * 4., settings.candidates});
        if (settings.candidates > 0)
            steps.push_back({settings.tolerance * 4., settings.candidates / 2});
    }
    DumpResult result;
    size_t     size = 0;
    for (const auto &[tolerance, count] : steps) {
        const Built built = build(objects, tolerance, count);
        std::string svg   = assemble(built, fingerprint, tolerance, count);
        size              = svg.size();
        if (settings.limit == 0 || size <= settings.limit) {
            result.svg        = std::move(svg);
            result.tolerance  = tolerance;
            result.candidates = count;
            result.dropped    = built.dropped;
            return result;
        }
    }
    result.error = "seamdump: " + std::to_string(size) + " bytes, above the limit of " + std::to_string(settings.limit) + " bytes";
    return result;
}

} // namespace PreciseSeam
} // namespace Slic3r
