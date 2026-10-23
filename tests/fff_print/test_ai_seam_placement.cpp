#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <optional>
#include <vector>

#include <boost/nowide/fstream.hpp>
#include <boost/filesystem/operations.hpp>

#include "test_helpers.hpp"
#include "test_utils.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/GCode/AISeamPlacement.hpp"
#include "libslic3r/GCode/AISeamWorkspace.hpp"
#include "libslic3r/Exception.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Utils.hpp"

using namespace Slic3r;
using namespace Slic3r::AISeamPlacement;

namespace {
struct WorkspaceDirectory {
    boost::filesystem::path path = boost::filesystem::temp_directory_path() /
        boost::filesystem::unique_path("orca-seam-test-%%%%-%%%%-%%%%");
    WorkspaceDirectory() { boost::filesystem::create_directory(path); }
    ~WorkspaceDirectory()
    {
        // Only remove the unique temporary directory created by this test instance.
        boost::system::error_code error;
        boost::filesystem::remove_all(path, error);
    }
};

Json read_json(const std::string &path)
{
    boost::nowide::ifstream input(path);
    return Json::parse(input);
}

// Exact integer export written beside every AI request; tests use it for precise geometry.
Json read_raw(const std::string &request_path)
{
    boost::filesystem::path raw(request_path);
    raw.replace_extension();
    raw += "-raw.json";
    return read_json(raw.string());
}

// An export target whose -raw.json sibling is removed together with it.
class ExportFile
{
public:
    boost::filesystem::path path() const { return m_dir.path() / "export.json"; }
    std::string string() const { return path().string(); }

private:
    ScopedTemporaryDir m_dir{"orca-seam-export"};
};

// The loop segment nearest to point, as the response "edge" pair of path points in millimeters.
Json edge_of(const Json &loop, double unit, const Vec2d &point)
{
    Json best;
    double best_distance = std::numeric_limits<double>::infinity();
    for (const Json &path : loop["paths"])
        for (size_t i = 0; i + 1 < path["points"].size(); ++i) {
            const Vec2d a(path["points"][i][0].get<double>() * unit, path["points"][i][1].get<double>() * unit);
            const Vec2d b(path["points"][i + 1][0].get<double>() * unit, path["points"][i + 1][1].get<double>() * unit);
            const Vec2d edge = b - a;
            const double t = edge.squaredNorm() == 0. ? 0. : std::clamp((point - a).dot(edge) / edge.squaredNorm(), 0., 1.);
            const double distance = (point - (a + t * edge)).norm();
            if (distance < best_distance) {
                best_distance = distance;
                best = Json::array({Json::array({a.x(), a.y()}), Json::array({b.x(), b.y()})});
            }
        }
    return best;
}

// Places one occurrence's seam at point, naming the loop segment it lies on.
void set_seam(Json &response, const Json &loop, size_t occurrence, double unit, const Vec2d &point)
{
    for (Json &seam : response["seams"])
        if (seam["occurrence_id"] == loop["occurrences"][occurrence]["occurrence_id"]) {
            seam["x_mm"] = point.x();
            seam["y_mm"] = point.y();
            seam["edge"] = edge_of(loop, unit, point);
        }
}

Json echoed_response(const Json &raw)
{
    Json result = {{"format", "orca-ai-seams-response"}, {"version", 1},
        {"export_id", raw["export_id"]}, {"export_digest", raw["export_digest"]},
        {"plate", raw["plate"]}, {"seams", Json::array()}};
    const double unit = raw["coord_unit_mm"];
    for (const Json &object : raw["objects"])
        for (const Json &layer : object["layers"])
            for (const Json &loop : layer["loops"])
                for (const Json &occurrence : loop["occurrences"]) {
                    const Vec2d seam(occurrence["source_seam"][0].get<double>() * unit,
                                     occurrence["source_seam"][1].get<double>() * unit);
                    result["seams"].push_back({{"occurrence_id", occurrence["occurrence_id"]},
                        {"x_mm", seam.x()}, {"y_mm", seam.y()}, {"edge", edge_of(loop, unit, seam)}});
                }
    return result;
}

void slice_cube(Print &print, Model &model)
{
    // Rear avoids expensive visibility sampling; two walls exercise the unchanged internal source.
    Test::init_print({Test::cube(2.)}, print, model,
        {{"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
         {"wall_loops", "2"}, {"seam_position", "back"}, {"skirt_loops", "0"}});
    print.process();
}
} // namespace

TEST_CASE("An isolated export preserves the connected response and working Gcode state", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    const bool done_before = print.is_step_done(psGCodeExport);
    ExportFile exported;
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    REQUIRE(print.is_step_done(psGCodeExport) == done_before);
    REQUIRE(print.ai_seam_options.response_path.empty());
    const Json output = read_raw(exported.string());
    REQUIRE(output["summary"]["occurrence_count"].get<size_t>() > 0);
    REQUIRE_FALSE(output.contains("environment"));
    REQUIRE(output["objects"][0]["layers"][0]["loops"][0].contains("precise_seam"));

    ScopedTemporaryFile response_file(".json");
    const Json response = echoed_response(output);
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    print.set_ai_seam_response(response_file.string());
    REQUIRE_NOTHROW(Test::gcode(print));
    REQUIRE(print.ai_seam_status == "applied");
    print.export_ai_seams(options);
    REQUIRE(print.ai_seam_options.response_path == response_file.string());
    REQUIRE(print.ai_seam_status == "applied");
}

TEST_CASE("The AI request uses rounded millimeters and its points are accepted", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json request = read_json(exported.string());
    const Json raw = read_raw(exported.string());
    REQUIRE_FALSE(request.contains("coord_unit_mm"));
    REQUIRE(request["precision_mm"].get<double>() == 0.01);
    REQUIRE(request["export_id"] == raw["export_id"]);
    REQUIRE(request["export_digest"] == raw["export_digest"]);
    const double unit = raw["coord_unit_mm"];
    const double raw_x = raw["objects"][0]["layers"][0]["loops"][0]["paths"][0]["points"][0][0].get<double>() * unit;
    const double x = request["objects"][0]["layers"][0]["loops"][0]["paths"][0]["points"][0][0].get<double>();
    REQUIRE(std::abs(x - raw_x) <= 0.005 + 1e-9);
    REQUIRE(std::abs(x * 100. - std::round(x * 100.)) < 1e-6);
    // Echo the rounded source seams exactly as an AI reading only the request would.
    Json response = {{"format", "orca-ai-seams-response"}, {"version", 1}, {"export_id", request["export_id"]},
        {"export_digest", request["export_digest"]}, {"plate", request["plate"]}, {"seams", Json::array()}};
    for (const Json &layer : request["objects"][0]["layers"])
        for (const Json &loop : layer["loops"])
            for (const Json &occurrence : loop["occurrences"]) {
                const Vec2d seam(occurrence["source_seam"][0].get<double>(), occurrence["source_seam"][1].get<double>());
                response["seams"].push_back({{"occurrence_id", occurrence["occurrence_id"]},
                    {"x_mm", seam.x()}, {"y_mm", seam.y()}, {"edge", edge_of(loop, 1., seam)}});
            }
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    print.set_ai_seam_response(response_file.string());
    REQUIRE_NOTHROW(Test::gcode(print));
    REQUIRE(print.ai_seam_status == "applied");
}

TEST_CASE("Workspace export is explicit and separates each plate's files", "[AISeamPlacement][AISeamWorkspace]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    WorkspaceDirectory folder;
    WorkspaceSettings settings;
    bool generated = false;
    const auto generate = [&] { generated = true; Test::gcode(print); };
    REQUIRE(generate_in_workspace(print, settings, Json::object(), generate).state == WorkspaceState::Idle);
    REQUIRE(generated);
    settings.enabled = true;
    generated = false;
    REQUIRE_THROWS_AS(generate_in_workspace(print, settings, Json::object(), generate), SlicingError);
    REQUIRE_FALSE(generated);
    settings.directory = folder.path.string();
    print.set_plate_index(0);
    print.set_ai_seam_response({});
    REQUIRE(generate_in_workspace(print, settings, {{"executable", "private-path"}}, generate).state == WorkspaceState::Exported);
    const auto first = folder.path / workspace_filename(0, false);
    REQUIRE(boost::filesystem::is_regular_file(first));
    REQUIRE(read_raw(first.string())["export_id"] == read_json(first.string())["export_id"]);
    REQUIRE_FALSE(read_json(first.string()).contains("environment"));
    settings.apply_response = false;
    print.set_ai_seam_response({});
    REQUIRE(generate_in_workspace(print, settings, Json::object(), generate).state == WorkspaceState::Exported);
    settings.apply_response = true;
    print.set_plate_index(1);
    print.set_ai_seam_response({});
    settings.share_paths = true;
    REQUIRE(generate_in_workspace(print, settings, {{"executable", "shared-path"}}, generate).state == WorkspaceState::Exported);
    const Json second = read_json((folder.path / workspace_filename(1, false)).string());
    REQUIRE(second["environment"]["executable"] == "shared-path");
    REQUIRE(boost::filesystem::is_regular_file(first));
}

TEST_CASE("Workspace snapshots responses and can recover from a corrupt answer", "[AISeamPlacement][AISeamWorkspace]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    print.set_plate_index(0);
    WorkspaceDirectory folder;
    WorkspaceSettings settings;
    settings.enabled = true;
    settings.directory = folder.path.string();
    const auto generate = [&] { Test::gcode(print); };
    REQUIRE(generate_in_workspace(print, settings, Json::object(), generate).state == WorkspaceState::Exported);
    const auto response_path = folder.path / workspace_filename(0, true);
    const Json response = echoed_response(read_raw((folder.path / workspace_filename(0, false)).string()));
    REQUIRE_FALSE(bool(write_file_atomically(response_path.string(), response.dump())));
    print.set_ai_seam_response({});
    // Replacing the public answer after preflight cannot change this job's frozen input.
    const auto mutate_original = [&] {
        REQUIRE_FALSE(bool(write_file_atomically(response_path.string(), "broken")));
        Test::gcode(print);
    };
    REQUIRE(generate_in_workspace(print, settings, Json::object(), mutate_original).state == WorkspaceState::Applied);
    REQUIRE(print.ai_seam_options.response_path.empty());
    print.set_ai_seam_response({});
    bool called = false;
    REQUIRE_THROWS_AS(generate_in_workspace(print, settings, Json::object(), [&] { called = true; }), SlicingError);
    REQUIRE_FALSE(called);
    settings.apply_response = false;
    REQUIRE(generate_in_workspace(print, settings, Json::object(), generate).state == WorkspaceState::Ignored);
    boost::nowide::ifstream original(response_path.string());
    std::string text;
    original >> text;
    REQUIRE(text == "broken");
}

TEST_CASE("Workspace geometry mismatch carries notes into a new request", "[AISeamPlacement][AISeamWorkspace]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    print.set_plate_index(0);
    WorkspaceDirectory folder;
    WorkspaceSettings settings;
    settings.enabled = true;
    settings.directory = folder.path.string();
    const auto generate = [&] { Test::gcode(print); };
    generate_in_workspace(print, settings, Json::object(), generate);
    const auto request_path = folder.path / workspace_filename(0, false);
    Json response = echoed_response(read_raw(request_path.string()));
    response["export_digest"] = std::string(64, '0');
    response["playbook"] = {{"format", "orca-ai-seams-playbook"}, {"version", 1}, {"text", "Keep seams on the rear."}};
    const auto response_path = folder.path / workspace_filename(0, true);
    REQUIRE_FALSE(bool(write_file_atomically(response_path.string(), response.dump())));
    print.set_ai_seam_response({});
    REQUIRE(generate_in_workspace(print, settings, Json::object(), generate).state == WorkspaceState::GeometryChanged);
    REQUIRE(read_json(request_path.string())["playbook"] == response["playbook"]);
    REQUIRE(read_json(response_path.string()) == response);
}

TEST_CASE("Matching geometry requires a complete near-boundary response", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    Json response = echoed_response(read_raw(exported.string()));
    response["seams"].erase(response["seams"].begin());
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    REQUIRE_THROWS(print.set_ai_seam_response(response_file.string()));
    response = echoed_response(read_raw(exported.string()));
    response["seams"][0]["x_mm"] = 10000.;
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    REQUIRE_THROWS(print.set_ai_seam_response(response_file.string()));
    REQUIRE(print.ai_seam_options.response_path.empty());
}

TEST_CASE("Layout coordinates in a response are named in the rejection", "[AISeamPlacement]")
{
    Print print;
    Model model;
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("seam_position", "back");
    config.set_deserialize_strict("layer_height", "0.2");
    config.set_deserialize_strict("initial_layer_print_height", "0.2");
    // Two arranged copies: at least one instance has a nonzero shift.
    Test::init_print(std::vector<TriangleMesh>{Test::cube(2.)}, print, model, config, nullptr, true, 2);
    print.process();
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json request = read_raw(exported.string());
    Json response = echoed_response(request);
    // Add each occurrence's own instance shift, as an AI returning layout XY would do.
    const double unit = request["coord_unit_mm"];
    std::map<std::string, Vec2d> shifts;
    for (const Json &instance : request["objects"][0]["instances"])
        shifts[instance["instance_key"]] = Vec2d(instance["shift"][0].get<double>(), instance["shift"][1].get<double>()) * unit;
    REQUIRE(std::any_of(shifts.begin(), shifts.end(), [](const auto &entry) { return entry.second.norm() > 0.1; }));
    std::map<std::string, std::string> instance_of;
    for (const Json &layer : request["objects"][0]["layers"])
        for (const Json &loop : layer["loops"])
            for (const Json &occurrence : loop["occurrences"])
                instance_of[occurrence["occurrence_id"]] = occurrence["instance_key"];
    for (Json &seam : response["seams"]) {
        const Vec2d &shift = shifts.at(instance_of.at(seam["occurrence_id"]));
        seam["x_mm"] = seam["x_mm"].get<double>() + shift.x();
        seam["y_mm"] = seam["y_mm"].get<double>() + shift.y();
        for (Json &end : seam["edge"])
            end = Json::array({end[0].get<double>() + shift.x(), end[1].get<double>() + shift.y()});
    }
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    REQUIRE_THROWS_WITH(print.set_ai_seam_response(response_file.string()),
        Catch::Matchers::ContainsSubstring("layout coordinates"));
    REQUIRE(print.ai_seam_options.response_path.empty());
}

TEST_CASE("The exact raw log never overwrites the connected response", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    ExportFile exported, request;
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    // The connected answer sits exactly where the next export's -raw.json would be published.
    const std::string response_path = raw_export_path(request.string());
    const std::string bytes = echoed_response(read_raw(exported.string())).dump();
    REQUIRE_FALSE(bool(write_file_atomically(response_path, bytes)));
    print.set_ai_seam_response(response_path);
    options.export_path = request.string();
    REQUIRE_THROWS_WITH(print.export_ai_seams(options), Catch::Matchers::ContainsSubstring("connected response"));
    REQUIRE_FALSE(boost::filesystem::exists(request.path()));
    REQUIRE(print.ai_seam_options.response_path == response_path);
    boost::nowide::ifstream input(response_path, std::ios::binary);
    REQUIRE(std::string(std::istreambuf_iterator<char>(input), {}) == bytes);
}

TEST_CASE("A response point is accepted within 0.05 mm of its exact loop", "[AISeamPlacement]")
{
    const double offset = GENERATE(0.04, 0.06);
    Print print;
    Model model;
    slice_cube(print, model);
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json raw = read_raw(exported.string());
    Json response = echoed_response(raw);
    // Move one seam perpendicular to the middle of a cube wall; other walls stay far away.
    const double unit = raw["coord_unit_mm"];
    const Json &loop = raw["objects"][0]["layers"][0]["loops"][0];
    const Json &points = loop["paths"][0]["points"];
    const Vec2d a(points[0][0].get<double>() * unit, points[0][1].get<double>() * unit);
    const Vec2d b(points[1][0].get<double>() * unit, points[1][1].get<double>() * unit);
    REQUIRE((b - a).norm() > 0.5);
    const Vec2d edge = (b - a).normalized();
    const Vec2d moved = (a + b) * 0.5 + Vec2d(-edge.y(), edge.x()) * offset;
    set_seam(response, loop, 0, unit, moved);
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    if (offset < 0.05)
        REQUIRE_NOTHROW(print.set_ai_seam_response(response_file.string()));
    else
        REQUIRE_THROWS_WITH(print.set_ai_seam_response(response_file.string()), Catch::Matchers::ContainsSubstring("too far"));
}

TEST_CASE("The response edge selects its segment and falls back visibly", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json raw = read_raw(exported.string());
    const double unit = raw["coord_unit_mm"];
    const Json &loop = raw["objects"][0]["layers"][0]["loops"][0];
    const Json &points = loop["paths"][0]["points"];
    REQUIRE(points.size() >= 4);
    const auto at = [&](size_t i) { return Json::array({points[i][0].get<double>() * unit, points[i][1].get<double>() * unit}); };
    const Vec2d middle((points[0][0].get<double>() + points[1][0].get<double>()) * 0.5 * unit,
                       (points[0][1].get<double>() + points[1][1].get<double>()) * 0.5 * unit);
    // Writes a response whose first seam sits mid-wall with the given edge; returns its fallback count.
    const auto fallbacks = [&](const Json &edge, bool strict = false) {
        Json response = echoed_response(raw);
        set_seam(response, loop, 0, unit, middle);
        for (Json &seam : response["seams"])
            if (seam["occurrence_id"] == loop["occurrences"][0]["occurrence_id"]) seam["edge"] = edge;
        REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
        Options applying;
        applying.response_path = response_file.string();
        applying.strict = strict;
        return Session(print, applying).edge_fallbacks().size();
    };
    REQUIRE(fallbacks(Json::array({at(0), at(1)})) == 0);
    REQUIRE(fallbacks(Json::array({at(1), at(0)})) == 0); // Swapped ends of the same segment.
    // Another wall's edge or an unknown end: the point is on its loop, so it is placed and reported.
    REQUIRE(fallbacks(Json::array({at(2), at(3)})) == 1);
    REQUIRE(fallbacks(Json::array({Json::array({100., 100.}), at(1)})) == 1);
    REQUIRE_THROWS_WITH(fallbacks(Json::array({at(2), at(3)}), true), Catch::Matchers::ContainsSubstring("too far from its edge"));
    REQUIRE_THROWS_WITH(fallbacks(Json::array({Json::array({100., 100.}), at(1)}), true),
        Catch::Matchers::ContainsSubstring("edge is not a segment"));
}

TEST_CASE("Geometry mismatch falls back for the whole plate and strict mode rejects it", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    ScopedTemporaryFile response_file(".json");
    const Json stale = {{"format", "orca-ai-seams-response"}, {"version", 1}, {"export_id", "stale"},
        {"export_digest", std::string(64, '0')}, {"seams", Json::array()}};
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), stale.dump())));
    Options options;
    options.response_path = response_file.string();
    Session baseline(print, options);
    REQUIRE_FALSE(baseline.applying());
    REQUIRE_FALSE(baseline.warning().empty());
    options.strict = true;
    REQUIRE_THROWS(Session(print, options));
}

TEST_CASE("Copies of one source loop can receive different finalized seams", "[AISeamPlacement]")
{
    Print print;
    Model model;
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("seam_position", "back");
    config.set_deserialize_strict("layer_height", "0.2");
    config.set_deserialize_strict("initial_layer_print_height", "0.2");
    config.set_deserialize_strict("wall_loops", "2");
    Test::init_print(std::vector<TriangleMesh>{Test::cube(2.)}, print, model, config, nullptr, true, 2);
    print.process();
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_raw(exported.string());
    Json response = echoed_response(output);
    const Json &loop = output["objects"][0]["layers"][0]["loops"][0];
    REQUIRE(loop["occurrences"].size() == 2);
    const double unit = output["coord_unit_mm"];
    for (size_t i = 0; i < 2; ++i) {
        const Json &point = loop["paths"][0]["points"][i];
        set_seam(response, loop, i, unit, Vec2d(point[0].get<double>() * unit, point[1].get<double>() * unit));
    }
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    print.set_ai_seam_response(response_file.string());
    Session applied(print, print.ai_seam_options);
    const auto *source = applied.sources()[0].loop;
    REQUIRE(applied.final_position(*source, 0)->point != applied.final_position(*source, 1)->point);
    REQUIRE_NOTHROW(Test::gcode(print));
    REQUIRE(print.ai_seam_status == "applied");
}

TEST_CASE("Different wall filaments preserve source-loop identity through clones", "[AISeamPlacement]")
{
    Print print;
    Model model;
    const auto config = Test::multifilament_config(2, {{"seam_position", "back"},
        {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"}, {"wall_loops", "2"},
        {"outer_wall_filament_id", "1"}, {"inner_wall_filament_id", "2"}, {"enable_prime_tower", "0"}});
    Test::init_print({Test::cube(2.)}, print, model, config);
    print.process();
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    REQUIRE_NOTHROW(print.export_ai_seams(options));
    const Json response = echoed_response(read_raw(exported.string()));
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    print.set_ai_seam_response(response_file.string());
    REQUIRE_NOTHROW(Test::gcode(print));
    REQUIRE(print.ai_seam_status == "applied");
}

TEST_CASE("Independent prints of the same model keep the same geometry digest", "[AISeamPlacement]")
{
    Print first, second;
    Model first_model, second_model;
    slice_cube(first, first_model);
    slice_cube(second, second_model);
    Options options;
    options.export_path = "unused-for-registry.json";
    Session a(first, options), b(second, options);
    REQUIRE(a.digest() == b.digest());
}

TEST_CASE("Playbook decisions survive new layer geometry without reusing old assignments", "[AISeamPlacement]")
{
    Print first, second;
    Model first_model, second_model;
    slice_cube(first, first_model);
    ExportFile exported, next_export;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    first.export_ai_seams(options);
    const Json old_export = read_raw(exported.string());
    Json response = echoed_response(old_export);
    const Json playbook = {{"format", "orca-ai-seams-playbook"}, {"version", 1},
        {"text", "Keep seams on the rear side. Recheck overhang policy after a layer-height change."}};
    response["playbook"] = playbook;
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    Test::init_print({Test::cube(2.)}, second, second_model,
        {{"layer_height", "0.25"}, {"initial_layer_print_height", "0.25"}, {"wall_loops", "2"}, {"seam_position", "back"}});
    second.process();
    second.set_ai_seam_response(response_file.string());
    options.export_path = next_export.string();
    second.export_ai_seams(options);
    const Json next = read_json(next_export.string());
    REQUIRE(next["playbook"] == playbook);
    REQUIRE(next["export_digest"] != old_export["export_digest"]);
    REQUIRE(second.ai_seam_options.response_path == response_file.string());
    options.export_path = response_file.string();
    REQUIRE_THROWS(second.export_ai_seams(options));
    REQUIRE(read_json(response_file.string())["format"] == "orca-ai-seams-response");
}

TEST_CASE("A late missing occurrence publishes neither Gcode nor a new seam export", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    ExportFile exported;
    ScopedTemporaryFile gcode_file(".gcode");
    Options options;
    options.export_path = exported.string();
    Session registry(print, options);
    REQUIRE_FALSE(registry.sources().empty());
    auto *source = const_cast<ExtrusionLoop *>(registry.sources()[0].loop);
    bool changed = false;
    size_t progress_events = 0;
    print.set_status_callback([&](const PrintBase::SlicingStatus &status) {
        // The serial G-code layer callback occurs after registry construction and before emission.
        // The first 80% event starts export; the next is the first generated layer, in any locale.
        if (status.percent == 80 && ++progress_events == 2) {
            source->set_loop_role(elrSkirt);
            changed = true;
        }
    });
    print.ai_seam_options = options;
    REQUIRE_THROWS(print.export_gcode(gcode_file.string(), nullptr));
    REQUIRE(changed);
    REQUIRE_FALSE(boost::filesystem::exists(gcode_file.path()));
    REQUIRE_FALSE(boost::filesystem::exists(gcode_file.string() + ".tmp"));
    REQUIRE_FALSE(boost::filesystem::exists(exported.path()));
    REQUIRE_FALSE(boost::filesystem::exists(raw_export_path(exported.string())));
}

TEST_CASE("Holes retain their own assignments in a complete export", "[AISeamPlacement]")
{
    Print print;
    Model model;
    Test::init_print({Test::TestMesh::cube_with_hole}, print, model,
        {{"seam_position", "back"}, {"wall_loops", "2"}, {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"}});
    print.process();
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_raw(exported.string());
    bool found_hole = false;
    for (const Json &object : output["objects"])
        for (const Json &layer : object["layers"])
            for (const Json &loop : layer["loops"]) found_hole = found_hole || loop["kind"] == "hole";
    REQUIRE(found_hole);
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), echoed_response(output).dump())));
    print.set_ai_seam_response(response_file.string());
    REQUIRE_NOTHROW(Test::gcode(print));
}

TEST_CASE("Fully overhanging external loops remain in the generation registry", "[AISeamPlacement]")
{
    const std::string generator = GENERATE("classic", "arachne");
    Print print;
    Model model;
    Test::init_print({Test::cube(2.)}, print, model,
        {{"seam_position", "back"}, {"wall_loops", "2"}, {"layer_height", "0.2"},
         {"initial_layer_print_height", "0.2"}, {"wall_generator", generator}});
    print.process();
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    Session before(print, options);
    // Reproduce the generator's role-loss condition while retaining actual outer-wall inset provenance.
    for (const auto &source : before.sources())
        for (auto &path : const_cast<ExtrusionLoop *>(source.loop)->paths) path.set_extrusion_role(erOverhangPerimeter);
    Session after(print, options);
    REQUIRE(after.sources().size() == before.sources().size());
    print.export_ai_seams(options);
    const Json output = read_raw(exported.string());
    REQUIRE(output["summary"]["source_loop_count"].get<size_t>() == before.sources().size());
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), echoed_response(output).dump())));
    print.set_ai_seam_response(response_file.string());
    REQUIRE_NOTHROW(Test::gcode(print));
}

TEST_CASE("Strong Precise Seam targets are distinct from recorded instance seams", "[AISeamPlacement]")
{
    Print print;
    Model model;
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("seam_position", "back");
    config.set_deserialize_strict("wall_loops", "2");
    config.set_deserialize_strict("layer_height", "0.2");
    config.set_deserialize_strict("initial_layer_print_height", "0.2");
    Test::init_print({Test::cube(2.)}, print, model, config);
    auto modifier_mesh = make_cube(0.5, 1., 2.);
    modifier_mesh.translate(0.75, -0.5, 0.);
    ModelVolume *modifier = model.objects.front()->add_volume(modifier_mesh);
    modifier->set_type(ModelVolumeType::PRECISE_SEAM_CENTER);
    modifier->name = "Front seam";
    print.apply(model, config);
    print.process();
    ExportFile exported;
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_raw(exported.string());
    size_t applied = 0;
    for (const Json &object : output["objects"])
        for (const Json &layer : object["layers"])
            for (const Json &loop : layer["loops"])
                if (loop["precise_seam"]["strong"]["status"] == "applied") {
                    ++applied;
                    const Json &strong = loop["precise_seam"]["strong"];
                    REQUIRE(strong["mode"] == "center");
                    REQUIRE(strong["modifier_name"] == "Front seam");
                    REQUIRE(strong["target_point"].size() == 2);
                    REQUIRE_FALSE(strong["target_on_overhang_path"].get<bool>());
                    REQUIRE(loop["precise_seam"]["weak"]["status"] == "not_evaluated_strong_applied");
                    for (const Json &occurrence : loop["occurrences"]) {
                        REQUIRE(occurrence["source_seam_matches_target"].is_boolean());
                        REQUIRE(occurrence["source_seam_target_distance_mm"].get<double>() >= 0.);
                    }
                }
    REQUIRE(applied > 0);
}

TEST_CASE("Generated outer-wall extrusion begins at the approved segment point", "[AISeamPlacement]")
{
    Print print;
    Model model;
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("seam_position", "back");
    config.set_deserialize_strict("layer_height", "0.2");
    config.set_deserialize_strict("initial_layer_print_height", "0.2");
    config.set_deserialize_strict("wall_loops", "2");
    config.set_deserialize_strict("seam_slope_type", "none");
    config.set_deserialize_strict("use_relative_e_distances", "0");
    // Use a wide, single-layer mesh so two perimeter walls and an interior split point fit.
    Test::init_print({make_cube(10., 10., 0.2)}, print, model, config);
    print.process();
    ExportFile exported;
    ScopedTemporaryFile response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_raw(exported.string());
    Json response = echoed_response(output);
    const Json &object = output["objects"][0];
    const Json &loop = object["layers"][0]["loops"][0];
    const Json &points = loop["paths"][0]["points"];
    const double unit = output["coord_unit_mm"];
    const Vec2d local((points[0][0].get<double>() + points[1][0].get<double>()) * 0.5 * unit,
                     (points[0][1].get<double>() + points[1][1].get<double>()) * 0.5 * unit);
    set_seam(response, loop, 0, unit, local);
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    print.set_ai_seam_response(response_file.string());
    const std::string gcode = Test::gcode(print);
    const Vec2d expected = local + Vec2d(object["instances"][0]["shift"][0].get<double>() * unit,
                                        object["instances"][0]["shift"][1].get<double>() * unit) - print.get_plate_origin().head<2>();
    const std::string role_tag = GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role);
    const std::string outer_tag = role_tag + ExtrusionEntity::role_to_string(erExternalPerimeter);
    std::optional<Vec2d> actual;
    bool outer = false;
    GCodeReader reader;
    reader.apply_config(config);
    reader.parse_buffer(gcode, [&](GCodeReader &current, const GCodeReader::GCodeLine &line) {
        if (line.comment().find(role_tag) == 0) outer = line.comment() == std::string_view(outer_tag);
        if (!actual && outer && line.extruding(current) && line.dist_XY(current) > 0.)
            actual = Vec2d(current.x(), current.y());
    });
    REQUIRE(actual.has_value());
    // Include normal 0.001 mm output quantization and the reader's float coordinates.
    REQUIRE_THAT(actual->x(), Catch::Matchers::WithinAbs(expected.x(), 0.0006));
    REQUIRE_THAT(actual->y(), Catch::Matchers::WithinAbs(expected.y(), 0.0006));
}
