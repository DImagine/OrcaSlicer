#ifndef slic3r_SeamDumpPrint_hpp_
#define slic3r_SeamDumpPrint_hpp_

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <boost/filesystem/path.hpp>
#include "PreciseSeamRules.hpp"
#include "PreciseSeamRulesApply.hpp"

// seamdump of a print: the request of the plate, a snapshot of its final outer loops and paint taken during the
// G-code export, and the files written from it where the project folder is known. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class ModelObject;
class Print;
namespace PreciseSeam {

struct SeamdumpFile {
    std::string object; // With `each`: the object's part of the file name, unique on the plate; empty otherwise.
    std::string svg;
    size_t      model_object = 0; // ObjectID of the model object with `each`.
};

// The files of one plate, built from a snapshot; writing them needs nothing from the print.
struct SeamdumpJob {
    SeamdumpRequest           request;
    std::vector<SeamdumpFile> files;
};

// The seamdump request of a plate from its model objects in list order: the first one of a Precise Seam volume
// that is not off. `conflicts`, when given, receives the volumes whose request differs.
std::optional<SeamdumpRequest> seamdump_request(const std::vector<const ModelObject *> &objects,
                                                std::vector<ModifierDiagnostic> *conflicts = nullptr);

// Builds the files that the plate's seamdump asks for, from the final perimeters, or null without a request.
// Conflicting requests, refused objects (spiral vase, belt printer, an open outer wall) and files above the size
// limit are reported in `warnings`. Without `each`, a refused object refuses the whole plate.
std::shared_ptr<const SeamdumpJob> build_seamdump(const Print &print, std::vector<std::string> &warnings,
                                                  const std::function<void()> &throw_if_canceled);

// File name parts of objects for `each`, from the whole object list: forbidden characters become "_", trailing
// dots and spaces go, reserved and empty names get a "_" before them; names equal without letter case are
// numbered "Name (1)", "Name (2)", ... in list order, and a name still taken gets the next free number.
std::vector<std::string> seamdump_object_names(const std::vector<std::string> &names);

// File names of the job in the project folder, in the order of its files: the requested name, or
// "<project>.seamdump.svg"; with `each`, "<base>-<object>.svg".
std::vector<std::string> seamdump_file_names(const SeamdumpJob &job, const std::string &project_name);

// File names that a request of a plate asks for, from its model objects in list order, whether its export is
// done or not; with `each`, one per object, in the order of `objects`.
std::vector<std::string> seamdump_planned_names(const SeamdumpRequest &request, const std::vector<const ModelObject *> &objects,
                                                const std::string &project_name);

// Writes one file through a temporary file that then replaces it, so the old file stays until the new one is
// complete. False on a failure, which leaves no temporary file.
bool write_seamdump_file(const boost::filesystem::path &path, const std::string &svg);

// The name without letter case, also beyond ASCII, as file systems that ignore it compare names.
std::string seamdump_fold(const std::string &name);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_SeamDumpPrint_hpp_
