// seamdump of a print: request, snapshot and files. Design: docs/HLSD/precise-seam.md
#include "SeamDumpPrint.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <locale>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <boost/filesystem/operations.hpp>
#include <boost/locale.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/nowide/fstream.hpp>
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/I18N.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/libslic3r.h"
#include "PreciseSeamRules.hpp"
#include "PreciseSeamRulesApply.hpp"
#include "SeamDump.hpp"
#include "SeamPlacer.hpp"

namespace Slic3r {
namespace PreciseSeam {

namespace {

// Outer loops of an entity tree in its order: closed loops with an outer wall path, or of inset 0, which a wholly
// overhanging outer loop keeps. An open outer wall cannot be dumped and sets `open`.
void collect_loops(const ExtrusionEntity &entity, std::vector<const ExtrusionLoop *> &out, bool &open)
{
    if (entity.is_collection()) {
        for (const ExtrusionEntity *child : static_cast<const ExtrusionEntityCollection &>(entity).entities)
            collect_loops(*child, out, open);
        return;
    }
    if (entity.is_loop()) {
        const auto &loop     = static_cast<const ExtrusionLoop &>(entity);
        bool        external = loop.inset_idx == 0;
        for (const ExtrusionPath &path : loop.paths)
            external = external || path.role() == erExternalPerimeter;
        if (external)
            out.push_back(&loop);
        return;
    }
    if (entity.role() == erExternalPerimeter)
        open = true;
}

// The model layers of one print object with their outer loops in plate coordinates of the first instance, and
// the paint along them. False for an open outer wall.
bool snapshot(const PrintObject &object, DumpObject &out, const std::function<void()> &throw_if_canceled)
{
    const SeamPaint paint(object);
    out.paint_source = paint.digest();
    const Point     shift = object.instances().front().shift_without_plate_offset();
    for (size_t i = 1; i < object.instances().size(); ++i)
        out.instances.push_back(unscale(Point(object.instances()[i].shift - object.instances().front().shift)));
    for (size_t li = 0; li < object.layers().size(); ++li) {
        throw_if_canceled();
        const Layer &layer = *object.layers()[li];
        DumpLayer    dump;
        dump.number  = li + 1;
        dump.slice_z = layer.slice_z;
        dump.print_z = layer.print_z;
        dump.height  = layer.height;
        for (const LayerRegion *region : layer.regions()) {
            std::vector<const ExtrusionLoop *> loops;
            bool                               open = false;
            for (const ExtrusionEntity *entity : region->perimeters.entities)
                collect_loops(*entity, loops, open);
            if (open)
                return false;
            const double reach = region->flow(frExternalPerimeter).width() / 2.;
            for (const ExtrusionLoop *loop : loops) {
                DumpLoop d;
                d.polygon = loop->polygon();
                d.hole    = (loop->loop_role() & elrHole) != 0;
                d.reach   = reach;
                // Paint is measured in object coordinates; arc lengths do not change with the shift.
                d.paint = paint.runs(d.polygon, layer.slice_z, reach);
                d.polygon.translate(shift);
                dump.loops.push_back(std::move(d));
            }
        }
        out.layers.push_back(std::move(dump));
    }
    return true;
}

// ASCII only, for the fallback of seamdump_fold().
std::string ascii_lower(std::string text)
{
    for (char &c : text)
        if (c >= 'A' && c <= 'Z')
            c = char(c - 'A' + 'a');
    return text;
}

} // namespace

std::string seamdump_fold(const std::string &name)
{
    // Lower case of the whole of Unicode needs a UTF-8 locale; without one, ASCII still folds.
    try {
        static const std::locale utf8 = boost::locale::generator().generate("en_US.UTF-8");
        return boost::locale::to_lower(name, utf8);
    } catch (...) {
        return ascii_lower(name);
    }
}

std::optional<SeamdumpRequest> seamdump_request(const std::vector<const ModelObject *> &objects,
                                                std::vector<ModifierDiagnostic> *conflicts)
{
    std::vector<const ModelVolume *> volumes;
    std::vector<LayerRules>          rules;
    // Only a rules field may hold a request, so other volumes are not parsed.
    for (const ModelObject *object : objects)
        for (const ModelVolume *volume : object->volumes)
            if (volume->is_precise_seam() && !volume->precise_seam_rules.empty()) {
                volumes.push_back(volume);
                rules.push_back(parse_layer_rules(volume->name, volume->precise_seam_rules));
            }
    std::vector<const LayerRules *> pointers;
    for (const LayerRules &r : rules)
        pointers.push_back(&r);
    const SeamdumpChoice choice = choose_seamdump(pointers);
    if (conflicts)
        for (size_t i : choice.conflicts)
            conflicts->push_back({volumes[i], {RuleDiagnosticCode::SeamdumpConflict, {true, 0}, rules[i].seamdump_element, "seamdump"}});
    return choice.request;
}

std::shared_ptr<const SeamdumpJob> build_seamdump(const Print &print, std::vector<std::string> &warnings,
                                                  const std::function<void()> &throw_if_canceled)
{
    // Print objects of each model object, in the order of the object list and then of their first instance.
    std::vector<std::pair<const ModelObject *, std::vector<const PrintObject *>>> groups;
    for (const ModelObject *model_object : print.model().objects) {
        std::vector<const PrintObject *> objects;
        for (const PrintObject *object : print.objects())
            if (object->model_object() == model_object && !object->instances().empty())
                objects.push_back(object);
        if (objects.empty())
            continue;
        const auto first = [model_object](const PrintObject *object) {
            const auto &list = model_object->instances;
            return std::find(list.begin(), list.end(), object->instances().front().model_instance) - list.begin();
        };
        std::stable_sort(objects.begin(), objects.end(), [&first](const PrintObject *a, const PrintObject *b) { return first(a) < first(b); });
        groups.emplace_back(model_object, std::move(objects));
    }
    std::vector<const ModelObject *> model_objects;
    for (const auto &group : groups)
        model_objects.push_back(group.first);
    std::vector<ModifierDiagnostic>      conflicts;
    const std::optional<SeamdumpRequest> request = seamdump_request(model_objects, &conflicts);
    for (std::string &warning : rule_warnings(conflicts))
        warnings.push_back(std::move(warning));
    if (!request)
        return nullptr;
    if (print.config().spiral_mode.value || print.config().belt_printer.value) {
        // TRN: Precise Seam seamdump. Do not translate "seamdump".
        warnings.push_back(_u8L("Precise Seam: seamdump does not support the spiral vase mode or belt printers, no file was written."));
        return nullptr;
    }
    DumpSettings settings;
    settings.tolerance         = double(request->tolerance) / double(nanometres_per_mm);
    settings.candidates        = request->candidates;
    settings.limit             = request->limit;
    settings.throw_if_canceled = throw_if_canceled;

    // File names of `each` come from the whole list, so a refused object keeps the numbers of the others.
    std::vector<std::string> names;
    for (const ModelObject *object : model_objects)
        names.push_back(object->name);
    const std::vector<std::string> parts = seamdump_object_names(names);

    auto job     = std::make_shared<SeamdumpJob>();
    job->request = *request;
    std::vector<DumpObject>  plate;
    std::vector<std::string> refused;
    for (size_t g = 0; g < groups.size(); ++g) {
        const auto &[model_object, objects] = groups[g];
        std::vector<DumpObject> dumped;
        bool                    open = false;
        for (const PrintObject *object : objects) {
            DumpObject d;
            d.name = model_object->name;
            if (!snapshot(*object, d, throw_if_canceled)) {
                open = true;
                break;
            }
            dumped.push_back(std::move(d));
        }
        if (open) {
            refused.push_back(model_object->name);
            continue;
        }
        if (!request->each) {
            for (DumpObject &d : dumped)
                plate.push_back(std::move(d));
            continue;
        }
        const DumpResult result = write_seamdump(dumped, settings);
        if (result.svg.empty())
            // TRN: Precise Seam seamdump. %1% is the object name, %2% and %3% numbers of bytes. Do not translate "seamdump".
            warnings.push_back((boost::format(_u8L("Precise Seam: seamdump of \"%1%\" was not written: %2% bytes are above the limit of %3% bytes.")) %
                                model_object->name % result.size % request->limit).str());
        else
            job->files.push_back({parts[g], result.svg});
    }
    for (const std::string &name : refused)
        // TRN: Precise Seam seamdump. %1% is the object name. Do not translate "seamdump".
        warnings.push_back((boost::format(_u8L("Precise Seam: seamdump does not support the open outer walls of \"%1%\".")) % name).str());
    if (!request->each) {
        // A plate file without some of its objects would look complete: none is written.
        if (!refused.empty())
            return nullptr;
        const DumpResult result = write_seamdump(plate, settings);
        if (result.svg.empty()) {
            // TRN: Precise Seam seamdump. %1% and %2% are numbers of bytes. Do not translate "seamdump".
            warnings.push_back((boost::format(_u8L("Precise Seam: seamdump was not written: %1% bytes are above the limit of %2% bytes.")) %
                                result.size % request->limit).str());
            return nullptr;
        }
        job->files.push_back({{}, result.svg});
    }
    if (job->files.empty())
        return nullptr;
    return job;
}

std::vector<std::string> seamdump_object_names(const std::vector<std::string> &names)
{
    std::vector<std::string> bases;
    for (const std::string &name : names) {
        std::string base;
        for (char c : name)
            base.push_back(static_cast<unsigned char>(c) < 0x20 || c == 0x7f || std::string("<>:\"/\\|?*").find(c) != std::string::npos ? '_' : c);
        while (!base.empty() && (base.back() == '.' || base.back() == ' '))
            base.pop_back();
        if (base.empty() || !portable_file_name(base + ".svg"))
            base = "_" + base;
        bases.push_back(base);
    }
    // Names equal without letter case are all numbered, in list order.
    std::map<std::string, size_t> count;
    for (const std::string &base : bases)
        ++count[seamdump_fold(base)];
    std::map<std::string, size_t> next;
    std::set<std::string>         used;
    std::vector<std::string>      out;
    for (const std::string &base : bases) {
        const std::string key      = seamdump_fold(base);
        const bool        numbered = count[key] > 1;
        size_t           &number   = next[key];
        std::string       name     = numbered ? base + " (" + std::to_string(++number) + ")" : base;
        // A name that an earlier object took, e.g. one called "Name (1)": the next free number.
        for (size_t k = numbered ? number : 0; used.count(seamdump_fold(name)) != 0;)
            name = base + " (" + std::to_string(++k) + ")";
        used.insert(seamdump_fold(name));
        out.push_back(name);
    }
    return out;
}

std::vector<std::string> seamdump_file_names(const SeamdumpJob &job, const std::string &project_name)
{
    std::vector<std::string> out;
    const std::string       &requested = job.request.file ? *job.request.file : std::string();
    // The requested name without its extension is the base of every object file.
    const std::string base = job.request.file ? requested.substr(0, requested.size() - 4) : project_name + ".seamdump";
    for (const SeamdumpFile &file : job.files)
        out.push_back(job.request.each ? base + "-" + file.object + ".svg" : (job.request.file ? requested : base + ".svg"));
    return out;
}

bool write_seamdump_file(const boost::filesystem::path &path, const std::string &svg)
{
    // Paths take UTF-8, as Orca sets up boost::filesystem at start.
    const boost::filesystem::path temp = path.parent_path() / (path.filename().string() + ".tmp");
    bool                          written = false;
    {
        boost::nowide::ofstream out(temp.string(), std::ios::binary | std::ios::trunc);
        if (out) {
            out.write(svg.data(), std::streamsize(svg.size()));
            out.close();
            written = !out.fail();
        }
    }
    boost::system::error_code error;
    if (written)
        boost::filesystem::rename(temp, path, error);
    if (!written || error) {
        boost::filesystem::remove(temp, error);
        return false;
    }
    BOOST_LOG_TRIVIAL(info) << "[PreciseSeamRules] seamdump: wrote " << path.string() << " bytes=" << svg.size();
    return true;
}

} // namespace PreciseSeam
} // namespace Slic3r
