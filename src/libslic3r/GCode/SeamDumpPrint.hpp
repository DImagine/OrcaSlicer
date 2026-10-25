#ifndef slic3r_SeamDumpPrint_hpp_
#define slic3r_SeamDumpPrint_hpp_

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <boost/filesystem/path.hpp>
#include "PreciseSeamRules.hpp"

// seamdump of a print: the request of the plate, a snapshot of its final outer loops and paint taken during the
// G-code export, and the files written from it where the project folder is known. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class Print;
namespace PreciseSeam {

struct SeamdumpFile {
    std::string object; // Model object name with `each`; empty for the file of the plate.
    std::string svg;
};

// The files of one plate, built from a snapshot; writing them needs nothing from the print.
struct SeamdumpJob {
    SeamdumpRequest           request;
    std::vector<SeamdumpFile> files;
};

// Builds the files that the plate's seamdump asks for, from the final perimeters, or null without a request.
// Conflicting requests, refused objects (spiral vase, belt printer, an open outer wall) and files above the size
// limit are reported in `warnings`. Without `each`, a refused object refuses the whole plate.
std::shared_ptr<const SeamdumpJob> build_seamdump(const Print &print, std::vector<std::string> &warnings,
                                                  const std::function<void()> &throw_if_canceled);

// File names of objects for `each`, without the extension: forbidden characters become "_", trailing dots and
// spaces go, reserved and empty names get a "_" before them; names equal without letter case are numbered
// "Name (1)", "Name (2)", ... in list order, and a name still taken gets the next free number.
std::vector<std::string> seamdump_object_names(const std::vector<std::string> &names);

// Writes the files of the job into the project folder, each through a temporary file that then replaces it.
// The plate file is the requested name, or "<project>.seamdump.svg"; with `each`, "<base>-<object>.svg".
// Returns one message per file that failed.
std::vector<std::string> write_seamdump_files(const SeamdumpJob &job, const boost::filesystem::path &folder,
                                              const std::string &project_name);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_SeamDumpPrint_hpp_
