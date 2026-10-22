#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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
Json read_json(const std::string &path)
{
    boost::nowide::ifstream input(path);
    return Json::parse(input);
}

Json echoed_response(const Json &exported)
{
    Json result = {{"format", "orca-ai-seams-response"}, {"version", 1},
        {"export_id", exported["export_id"]}, {"export_digest", exported["export_digest"]},
        {"coord_unit_mm", exported["coord_unit_mm"]}, {"plate", exported["plate"]}, {"seams", Json::array()}};
    const double unit = exported["coord_unit_mm"];
    for (const Json &object : exported["objects"])
        for (const Json &layer : object["layers"])
            for (const Json &loop : layer["loops"])
                for (const Json &occurrence : loop["occurrences"])
                    result["seams"].push_back({{"occurrence_id", occurrence["occurrence_id"]},
                        {"x_mm", occurrence["source_seam"][0].get<double>() * unit},
                        {"y_mm", occurrence["source_seam"][1].get<double>() * unit}});
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
    ScopedTemporaryFile exported(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    REQUIRE(print.is_step_done(psGCodeExport) == done_before);
    REQUIRE(print.ai_seam_options.response_path.empty());
    const Json output = read_json(exported.string());
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

TEST_CASE("Matching geometry requires a complete near-boundary response", "[AISeamPlacement]")
{
    Print print;
    Model model;
    slice_cube(print, model);
    ScopedTemporaryFile exported(".json"), response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    Json response = echoed_response(read_json(exported.string()));
    response["seams"].erase(response["seams"].begin());
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    REQUIRE_THROWS(print.set_ai_seam_response(response_file.string()));
    response = echoed_response(read_json(exported.string()));
    response["seams"][0]["x_mm"] = 10000.;
    REQUIRE_FALSE(bool(write_file_atomically(response_file.string(), response.dump())));
    REQUIRE_THROWS(print.set_ai_seam_response(response_file.string()));
    REQUIRE(print.ai_seam_options.response_path.empty());
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
    ScopedTemporaryFile exported(".json"), response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_json(exported.string());
    Json response = echoed_response(output);
    const Json &loop = output["objects"][0]["layers"][0]["loops"][0];
    REQUIRE(loop["occurrences"].size() == 2);
    const double unit = output["coord_unit_mm"];
    for (size_t i = 0; i < 2; ++i) {
        const Json &point = loop["paths"][0]["points"][i];
        for (Json &seam : response["seams"])
            if (seam["occurrence_id"] == loop["occurrences"][i]["occurrence_id"]) {
                seam["x_mm"] = point[0].get<double>() * unit;
                seam["y_mm"] = point[1].get<double>() * unit;
            }
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
    ScopedTemporaryFile exported(".json"), response_file(".json");
    Options options;
    options.export_path = exported.string();
    REQUIRE_NOTHROW(print.export_ai_seams(options));
    const Json response = echoed_response(read_json(exported.string()));
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
    ScopedTemporaryFile exported(".json"), response_file(".json"), next_export(".json");
    Options options;
    options.export_path = exported.string();
    first.export_ai_seams(options);
    const Json old_export = read_json(exported.string());
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
    ScopedTemporaryFile exported(".json"), gcode_file(".gcode");
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
}

TEST_CASE("Holes retain their own assignments in a complete export", "[AISeamPlacement]")
{
    Print print;
    Model model;
    Test::init_print({Test::TestMesh::cube_with_hole}, print, model,
        {{"seam_position", "back"}, {"wall_loops", "2"}, {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"}});
    print.process();
    ScopedTemporaryFile exported(".json"), response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_json(exported.string());
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
    ScopedTemporaryFile exported(".json"), response_file(".json");
    Options options;
    options.export_path = exported.string();
    Session before(print, options);
    // Reproduce the generator's role-loss condition while retaining actual outer-wall inset provenance.
    for (const auto &source : before.sources())
        for (auto &path : const_cast<ExtrusionLoop *>(source.loop)->paths) path.set_extrusion_role(erOverhangPerimeter);
    Session after(print, options);
    REQUIRE(after.sources().size() == before.sources().size());
    print.export_ai_seams(options);
    const Json output = read_json(exported.string());
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
    ScopedTemporaryFile exported(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_json(exported.string());
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
    ScopedTemporaryFile exported(".json"), response_file(".json");
    Options options;
    options.export_path = exported.string();
    print.export_ai_seams(options);
    const Json output = read_json(exported.string());
    Json response = echoed_response(output);
    const Json &object = output["objects"][0];
    const Json &loop = object["layers"][0]["loops"][0];
    const Json &points = loop["paths"][0]["points"];
    const double unit = output["coord_unit_mm"];
    const Vec2d local((points[0][0].get<double>() + points[1][0].get<double>()) * 0.5 * unit,
                     (points[0][1].get<double>() + points[1][1].get<double>()) * 0.5 * unit);
    for (Json &seam : response["seams"])
        if (seam["occurrence_id"] == loop["occurrences"][0]["occurrence_id"]) {
            seam["x_mm"] = local.x();
            seam["y_mm"] = local.y();
        }
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
