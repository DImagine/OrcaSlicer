#ifndef slic3r_SeamDump_hpp_
#define slic3r_SeamDump_hpp_

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"

// seamdump: the outer loops of every layer written to one SVG file for an external program, such as a language
// model, to choose the seams. Loops are simplified with a checked tolerance, numbered across layers and given
// seam candidates and the user's seam paint. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
namespace PreciseSeam {

// Seam paint along a source loop: arc lengths in mm from its first point in its direction; `end` may pass the
// loop length to wrap around its start.
struct DumpPaintRun {
    double begin   = 0.;
    double end     = 0.;
    bool   blocked = false; // Blocked paint; enforced otherwise. Blocked wins where both lie.
};

// One outer loop of a layer: the outer contour or a hole of an island.
struct DumpLoop {
    Polygon                   polygon; // Plate coordinates of the first instance, scaled.
    bool                      hole  = false;
    double                    reach = 0.; // Paint projection distance in mm, which the fingerprint records.
    std::vector<DumpPaintRun> paint;
};

// One model layer of a print object, numbered from 1; its loops in source order, which data-source-loop gives.
struct DumpLayer {
    size_t                number  = 0;
    double                slice_z = 0.;
    double                print_z = 0.;
    double                height  = 0.;
    std::vector<DumpLoop> loops;
};

struct DumpObject {
    std::string            name;
    std::vector<Vec2d>     instances;    // Shifts in mm of the other instances from the first one.
    std::string            paint_source; // Digest of the painted facets and their transformations.
    std::vector<DumpLayer> layers;       // Every model layer in order, with loops or without.
};

struct DumpSettings {
    double   tolerance  = 0.1; // T in mm.
    int      candidates = 8;   // K.
    uint64_t limit      = 0;   // Bytes; 0 is none.
    // Called between parts of the work, at least once per loop and layer; it throws to cancel.
    std::function<void()> throw_if_canceled;
};

struct DumpResult {
    std::string svg;      // Empty when no file fits the limit.
    size_t      size = 0; // Bytes of the last, coarsest file built when none fits.
    double      tolerance  = 0.; // T and K of the file written.
    int         candidates = 0;
    size_t      dropped    = 0; // Loops left out.
};

// The SVG of the objects, numbered O1, O2, ... in the given order. Above the limit, the file is built again
// with 2T, 4T and then half of K; if none fits, `svg` is empty.
DumpResult write_seamdump(const std::vector<DumpObject> &objects, const DumpSettings &settings);

// SHA-256 of the bytes as 64 hexadecimal digits.
std::string seamdump_sha256(const std::string &data);

// A number of hundredths of a millimetre as the file writes it: no trailing zeros, no exponent.
std::string seamdump_number(int64_t hundredths);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_SeamDump_hpp_
