#pragma once

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "libslic3r/Point.hpp"

namespace Slic3r {
class ExtrusionLoop;
class Layer;
class LayerRegion;
class Print;
class PrintObject;

namespace AISeamPlacement {
using Json = nlohmann::ordered_json;

// Runtime only: these options must never enter PrintConfig or the 3MF serializer.
struct Options {
    std::string response_path;
    std::string export_path;
    std::string prompt;
    Json environment = Json::object();
    Json playbook = nullptr; // Portable user decisions, independent of this geometry's IDs/points.
    bool strict = false;
    bool enabled() const { return !response_path.empty() || !export_path.empty(); }
};

struct Projection {
    size_t path_idx = 0;
    size_t segment_idx = 0;
    Point point = Point(0, 0);
    double parameter = 0.;
    double distance_mm = 0.;
    double s_mm = 0.;
    float unsupported_distance = 0.f; // Filled by preflight at the rounded point, zero on the first layer.
};

bool is_external(const ExtrusionLoop &loop);
Projection project(const ExtrusionLoop &loop, const Vec2d &point_mm);
// Uses the preflight segment/visit, bypassing all standard snapping and overhang preferences.
void split_exact(ExtrusionLoop &loop, const Projection &projection);
Json read_response(const std::string &path);
void validate_response(const Json &response);
std::string initial_prompt();

class Session {
public:
    struct Source {
        const PrintObject *object = nullptr;
        const Layer *layer = nullptr;
        const LayerRegion *region = nullptr;
        const ExtrusionLoop *loop = nullptr;
        std::string id;
        size_t object_idx = 0;
        size_t layer_idx = 0;
        Json annotations = Json::object();
    };

    Session(const Print &print, const Options &options);
    bool exporting() const { return !m_options.export_path.empty(); }
    bool applying() const { return m_applying; }
    const std::string &warning() const { return m_warning; }
    const std::string &digest() const { return m_digest; }
    const std::vector<Source> &sources() const { return m_sources; }
    void annotate(const ExtrusionLoop *loop, Json annotations);
    void add_alignment_group(Json group);
    void clear_clones() { m_clones.clear(); }
    void register_clone(const ExtrusionLoop &clone, const ExtrusionLoop &source);
    std::optional<std::string> source_id(const ExtrusionLoop *loop) const;
    std::optional<Projection> final_position(const ExtrusionLoop &source, size_t instance_idx) const;
    void record(const ExtrusionLoop &source, size_t instance_idx, const Point &actual_seam);
    void finish(bool write_export = false);

private:
    struct Occurrence {
        size_t source_idx = 0;
        size_t instance_idx = 0;
        std::string id;
        std::optional<Projection> target;
        std::optional<Point> seam;
    };
    const Print &m_print;
    Options m_options;
    std::vector<Source> m_sources;
    std::vector<Occurrence> m_occurrences;
    std::map<const ExtrusionLoop *, size_t> m_source_index;
    std::map<const ExtrusionLoop *, const ExtrusionLoop *> m_clones;
    std::map<std::pair<const ExtrusionLoop *, size_t>, size_t> m_occurrence_index;
    Json m_objects = Json::array();
    Json m_object_digests = Json::object();
    Json m_alignment_groups = Json::array();
    std::string m_digest;
    std::string m_warning;
    bool m_applying = false;
    const ExtrusionLoop &original(const ExtrusionLoop &loop) const;
};
} // namespace AISeamPlacement
} // namespace Slic3r
