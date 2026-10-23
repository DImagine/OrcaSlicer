#include "AISeamWorkspace.hpp"

#include <iterator>
#include <stdexcept>
#include <utility>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/nowide/cstdio.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/system/error_code.hpp>

#include "libslic3r/Exception.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintBase.hpp"

namespace Slic3r::AISeamPlacement {
namespace {
struct RestoreWorkspace {
    Print &print;
    Options previous;
    std::string temporary;
    ~RestoreWorkspace()
    {
        // The temporary response belongs to exactly one job, including failed/cancelled jobs.
        print.ai_seam_options = std::move(previous);
        if (!temporary.empty()) boost::nowide::remove(temporary.c_str());
    }
};
} // namespace

std::string workspace_filename(int plate_index, bool response)
{
    return "ai-seam-plate" + std::to_string(plate_index + 1) + (response ? "-response.json" : "-request.json");
}

WorkspaceResult generate_in_workspace(Print &print, const WorkspaceSettings &settings,
                                     const Json &environment, const std::function<void()> &generate)
{
    if (!settings.enabled) {
        generate();
        return {};
    }
    RestoreWorkspace restore{print, print.ai_seam_options, {}};
    print.ai_seam_status.clear();
    try {
        if (print.config().spiral_mode.value)
            throw std::runtime_error("AI Seam is unavailable in spiral vase mode");
        if (settings.directory.empty()) throw std::runtime_error("Choose an AI Seam working folder");
        const boost::filesystem::path folder = boost::filesystem::absolute(settings.directory);
        if (!boost::filesystem::is_directory(folder))
            throw std::runtime_error("AI Seam working folder does not exist: " + folder.string());
        const auto request_path = folder / workspace_filename(print.get_plate_index(), false);
        const auto response_path = folder / workspace_filename(print.get_plate_index(), true);
        boost::system::error_code equivalent_error;
        if (boost::filesystem::equivalent(request_path, response_path, equivalent_error))
            throw std::runtime_error("AI Seam request and response must be different files");

        // The UI has already invalidated G-code, so this cannot cancel its own running worker.
        print.set_ai_seam_response({});
        print.ai_seam_options = {};
        WorkspaceResult result;
        bool mismatch = false;
        boost::system::error_code existence_error;
        const bool response_exists = boost::filesystem::exists(response_path, existence_error);
        if (settings.apply_response && existence_error &&
            existence_error != boost::system::errc::no_such_file_or_directory)
            throw std::runtime_error("Cannot inspect AI Seam response: " + response_path.string());
        const bool has_response = settings.apply_response && response_exists;
        if (has_response) {
            boost::nowide::ifstream input(response_path.string(), std::ios::binary);
            if (!input) throw std::runtime_error("Cannot read AI Seam response: " + response_path.string());
            const std::string bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
            if (input.bad()) throw std::runtime_error("Cannot finish reading AI Seam response");
            restore.temporary = (boost::filesystem::temp_directory_path() /
                boost::filesystem::unique_path("orca-ai-response-%%%%-%%%%-%%%%.json")).string();
            {
                boost::nowide::ofstream output(restore.temporary, std::ios::binary | std::ios::trunc);
                output.write(bytes.data(), std::streamsize(bytes.size()));
                output.close();
                if (!output) throw std::runtime_error("Cannot snapshot AI Seam response");
            }
            const Json response = read_response(restore.temporary);
            print.set_ai_seam_response(restore.temporary);
            Session preflight(print, print.ai_seam_options);
            mismatch = !preflight.warning().empty();
            result.seam_count = response.at("seams").size();
        }
        if (!has_response || mismatch) {
            Options options;
            options.export_path = request_path.string();
            options.prompt = initial_prompt() + "\nSave the complete response as " +
                response_path.filename().string() + ". Do not overwrite the request file.\n";
            if (settings.share_paths) {
                options.environment = environment;
                options.environment["export_file"] = request_path.string();
                options.environment["response_file"] = response_path.string();
            }
            // Only a connected valid response contributes a playbook; ignored corrupt files stay unread.
            print.export_ai_seams(options);
        }
        generate();
        if (has_response && !mismatch && print.ai_seam_status != "applied")
            throw std::runtime_error("AI Seam response was not applied by G-code generation");
        result.state = mismatch ? WorkspaceState::GeometryChanged :
            has_response ? WorkspaceState::Applied :
            !settings.apply_response && response_exists ? WorkspaceState::Ignored : WorkspaceState::Exported;
        return result;
    } catch (const CanceledException &) {
        print.ai_seam_status.clear();
        throw;
    } catch (const SlicingError &) {
        print.ai_seam_status.clear();
        throw;
    } catch (const std::exception &error) {
        print.ai_seam_status.clear();
        // Ordinary exchange errors belong in the slicer's error notification, not its crash dialog.
        throw SlicingError(std::string("AI Seam: ") + error.what());
    }
}
} // namespace Slic3r::AISeamPlacement
