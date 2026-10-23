#pragma once

#include <cstddef>
#include <functional>
#include <string>

#include "AISeamPlacement.hpp"

namespace Slic3r { class Print; }
namespace Slic3r::AISeamPlacement {

// Session-only GUI settings: deliberately absent from profiles and project serialization.
struct WorkspaceSettings {
    bool enabled = false;
    bool apply_response = true;
    bool share_paths = false;
    std::string directory;
};

enum class WorkspaceState { Idle, Checking, Exported, Ignored, Applied, GeometryChanged, Error, Cancelled };
struct WorkspaceResult {
    WorkspaceState state = WorkspaceState::Idle;
    size_t seam_count = 0;
    size_t edge_fallback_count = 0; // Seams placed at the nearest loop point without a matching edge.
    std::string error;
};

std::string workspace_filename(int plate_index, bool response);
// Caller exclusively owns the sliced Print and invalidates G-code before starting its worker.
WorkspaceResult generate_in_workspace(Print &print, const WorkspaceSettings &settings,
                                     const Json &environment, const std::function<void()> &generate);
} // namespace Slic3r::AISeamPlacement
