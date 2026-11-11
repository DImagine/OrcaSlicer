#include <catch2/catch_all.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "libslic3r/Config.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/OuterWallPasses.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "test_helpers.hpp"

using namespace Slic3r;
using namespace Slic3r::Test;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace {

// Half of the layer height: the outer wall is split into two passes.
const double layer_height = 0.2;
const double pass_height  = 0.1;
// Middle of the allowed pass width range set below.
const double nominal_width = 0.5;

// Every setting the assertions depend on, so none of them rests on a default.
DynamicPrintConfig passes_config(const char *fallback)
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "split_outer_wall",            true },
        { "split_outer_wall_count",      2 },
        { "split_outer_wall_min_width",  "0.4" },
        { "split_outer_wall_max_width",  "0.6" },
        { "split_outer_wall_fallback",   fallback },
        { "split_outer_wall_holes",      true },
        { "wall_generator",              "classic" },
        { "nozzle_diameter",             "0.4" },
        { "layer_height",                layer_height },
        { "initial_layer_print_height",  layer_height },
        { "wall_loops",                  2 },
        { "sparse_infill_density",       "15%" },
        { "elefant_foot_compensation",   0. },
        { "skirt_loops",                 0 },
    });
    return config;
}

// Widths and heights of the paths of the passes of a layer, per pass, bottom pass first. The top pass is last.
struct PassPaths
{
    std::vector<double> widths;
    std::vector<double> heights;
};

void collect_passes(const ExtrusionEntityCollection &collection, PassPaths &out)
{
    for (const ExtrusionEntity *entity : collection.entities)
        if (entity->is_collection())
            collect_passes(*static_cast<const ExtrusionEntityCollection*>(entity), out);
        else if (entity->is_loop() && static_cast<const ExtrusionLoop*>(entity)->outer_wall_pass)
            for (const ExtrusionPath &path : static_cast<const ExtrusionLoop*>(entity)->paths) {
                out.widths.emplace_back(path.width);
                out.heights.emplace_back(path.height);
            }
}

// Number of the loops of passes in a collection.
size_t count_pass_loops(const ExtrusionEntityCollection &collection)
{
    size_t count = 0;
    for (const ExtrusionEntity *entity : collection.entities)
        if (entity->is_collection())
            count += count_pass_loops(*static_cast<const ExtrusionEntityCollection*>(entity));
        else if (entity->is_loop() && static_cast<const ExtrusionLoop*>(entity)->outer_wall_pass)
            ++ count;
    return count;
}

// What the G-code extrudes at which height.
struct ExtrusionLevels
{
    // Start point of the first extrusion at each height, by the height in microns.
    std::map<int, Vec2d> first_start;
    // Whether the nozzle ever extruded below a height it had already printed at.
    bool                 went_below { false };

    bool has(int z_microns) const { return this->first_start.count(z_microns) == 1; }
};

ExtrusionLevels extrusion_levels(const std::string &gcode)
{
    ExtrusionLevels   out;
    double            highest_z = 0.;
    GCodeReader       reader;
    reader.parse_buffer(gcode, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        if (line.extruding(self) && line.dist_XY(self) > 0.) {
            const double z = line.new_Z(self);
            // The reader still holds the position the move starts at.
            out.first_start.emplace(int(std::lround(z * 1000.)), Vec2d(self.x(), self.y()));
            if (z < highest_z - 1e-3)
                out.went_below = true;
            if (z > highest_z)
                highest_z = z;
        }
    });
    return out;
}

ExtrusionLevels extrusion_levels(Print &print)
{
    return extrusion_levels(Slic3r::Test::gcode(print));
}

std::vector<PassPaths> passes_of(const Layer &layer)
{
    std::vector<PassPaths> out;
    for (const LayerRegion *region : layer.regions()) {
        if (out.size() < region->outer_wall_lower_passes.size() + 1)
            out.resize(region->outer_wall_lower_passes.size() + 1);
        for (size_t pass = 0; pass < region->outer_wall_lower_passes.size(); ++ pass)
            collect_passes(region->outer_wall_lower_passes[pass], out[pass]);
        collect_passes(region->perimeters, out.back());
    }
    return out;
}

} // namespace

TEST_CASE("Outer wall of a vertical wall is printed in passes of the nominal width", "[OuterWallPasses]")
{
    Print print;
    init_and_process_print({ cube(20.) }, print, passes_config("split_layer"));
    const PrintObject &object = *print.objects().front();

    // Vertical walls fit any width range, so no layer is split: 20 mm at 0.2 mm.
    REQUIRE(object.layers().size() == 100);

    SECTION("the first layer keeps its regular outer wall") {
        const Layer &first = *object.layers().front();
        REQUIRE_FALSE(first.outer_wall_pass_slices);
        for (const PassPaths &pass : passes_of(first))
            REQUIRE(pass.widths.empty());
    }

    SECTION("the layers above print the outer wall in two passes") {
        for (size_t layer_idx : { size_t(1), size_t(50), size_t(99) }) {
            CAPTURE(layer_idx);
            const Layer &layer = *object.layers()[layer_idx];
            REQUIRE(layer.outer_wall_pass_slices);
            REQUIRE(layer.outer_wall_pass_slices->count() == 2);
            const std::vector<PassPaths> passes = passes_of(layer);
            REQUIRE(passes.size() == 2);
            for (const PassPaths &pass : passes) {
                REQUIRE_FALSE(pass.widths.empty());
                for (double width : pass.widths)
                    REQUIRE_THAT(width, WithinAbs(nominal_width, 0.005));
                for (double height : pass.heights)
                    REQUIRE_THAT(height, WithinAbs(pass_height, 1e-4));
            }
        }
    }
}

TEST_CASE("Passes of a sloped wall differ in width by the slope", "[OuterWallPasses]")
{
    // A cone of radius 10 and height 20 narrows by 0.5 mm per mm of height, that is by 0.05 mm per pass.
    Print print;
    init_and_process_print({ make_cone(10., 20.) }, print, passes_config("split_layer"));
    const PrintObject &object = *print.objects().front();

    const Layer &layer = *object.layers()[30];
    REQUIRE_THAT(layer.height, WithinAbs(layer_height, 1e-6));
    REQUIRE(layer.outer_wall_pass_slices);
    const std::vector<PassPaths> passes = passes_of(layer);
    REQUIRE(passes.size() == 2);
    REQUIRE_FALSE(passes.front().widths.empty());
    REQUIRE_FALSE(passes.back().widths.empty());
    // The bottom pass reaches further out than the top one.
    for (double width : passes.front().widths)
        REQUIRE_THAT(width, WithinAbs(nominal_width + 0.025, 0.01));
    for (double width : passes.back().widths)
        REQUIRE_THAT(width, WithinAbs(nominal_width - 0.025, 0.01));
}

TEST_CASE("Layers of a slope out of the pass width range follow the fallback setting", "[OuterWallPasses]")
{
    // A cone of radius 10 and height 2 narrows by 5 mm per mm of height, that is by 0.5 mm per pass:
    // more than the 0.2 mm the widths of the passes may differ by.
    SECTION("split layer: the layers are sliced at the pass height") {
        Print print;
        init_and_process_print({ make_cone(10., 2.) }, print, passes_config("split_layer"));
        const PrintObject &object = *print.objects().front();
        REQUIRE(object.layers().size() > 10);
        for (size_t layer_idx = 1; layer_idx < object.layers().size(); ++ layer_idx) {
            CAPTURE(layer_idx);
            const Layer &layer = *object.layers()[layer_idx];
            REQUIRE_FALSE(layer.outer_wall_pass_slices);
            REQUIRE_THAT(layer.height, WithinAbs(pass_height, 1e-6));
        }
    }
    SECTION("regular: the layers keep their height and their regular outer wall") {
        Print print;
        init_and_process_print({ make_cone(10., 2.) }, print, passes_config("regular"));
        const PrintObject &object = *print.objects().front();
        REQUIRE(object.layers().size() <= 10);
        for (size_t layer_idx = 0; layer_idx < object.layers().size(); ++ layer_idx) {
            CAPTURE(layer_idx);
            const Layer &layer = *object.layers()[layer_idx];
            REQUIRE_FALSE(layer.outer_wall_pass_slices);
            REQUIRE_THAT(layer.height, WithinAbs(layer_height, 1e-6));
            for (const PassPaths &pass : passes_of(layer))
                REQUIRE(pass.widths.empty());
        }
    }
}

TEST_CASE("Passes are printed below the layer before the rest of the layer", "[OuterWallPasses]")
{
    Print print;
    init_and_process_print({ cube(20.) }, print, passes_config("split_layer"));
    const ExtrusionLevels levels = extrusion_levels(print);

    // The nozzle never returns below a height it already printed at.
    REQUIRE_FALSE(levels.went_below);
    // The first layer is printed at its height only, the layers above also half a layer below.
    REQUIRE_FALSE(levels.has(100));
    REQUIRE(levels.has(200));
    REQUIRE(levels.has(300));
    REQUIRE(levels.has(400));
    REQUIRE(levels.has(19900));
    REQUIRE(levels.has(20000));
}

TEST_CASE("Passes of all the objects are printed lowest first", "[OuterWallPasses]")
{
    // The low passes of both cubes have to be down before the rest of the layer of either of them.
    Print print;
    init_and_process_print({ cube(20.), cube(20.) }, print, passes_config("split_layer"));
    const ExtrusionLevels levels = extrusion_levels(print);
    REQUIRE(levels.has(300));
    REQUIRE_FALSE(levels.went_below);
}

TEST_CASE("Odd number of passes takes the middle pass for the middle of the layer", "[OuterWallPasses]")
{
    // Three passes of 0.08 mm. The cone narrows by 0.5 mm per mm of height, that is by 0.04 mm per pass.
    DynamicPrintConfig config = passes_config("split_layer");
    config.set_deserialize_strict({
        { "split_outer_wall_count",      3 },
        { "layer_height",                0.24 },
        { "initial_layer_print_height",  0.24 },
    });
    Print print;
    init_and_process_print({ make_cone(10., 20.) }, print, config);
    const Layer &layer = *print.objects().front()->layers()[30];
    REQUIRE(layer.outer_wall_pass_slices);
    REQUIRE(layer.outer_wall_pass_slices->count() == 3);
    const std::vector<PassPaths> passes = passes_of(layer);
    REQUIRE(passes.size() == 3);
    const double expected[] = { nominal_width + 0.04, nominal_width, nominal_width - 0.04 };
    for (size_t pass = 0; pass < 3; ++ pass) {
        CAPTURE(pass);
        REQUIRE_FALSE(passes[pass].widths.empty());
        for (double width : passes[pass].widths)
            REQUIRE_THAT(width, WithinAbs(expected[pass], 0.01));
        for (double height : passes[pass].heights)
            REQUIRE_THAT(height, WithinAbs(0.08, 1e-4));
    }
}

TEST_CASE("All the passes of an outer wall start at the same point", "[OuterWallPasses]")
{
    // The seam placer decides a random seam per wall ahead of the G-code, so this guards the property itself,
    // whichever way the passes come to share the point.
    DynamicPrintConfig config = passes_config("split_layer");
    config.set_deserialize_strict({ { "seam_position", "random" } });
    Print print;
    init_and_process_print({ make_cone(10., 20.) }, print, config);
    const ExtrusionLevels levels = extrusion_levels(print);
    for (int layer_z : { 2000, 6200, 10400 }) {
        CAPTURE(layer_z);
        REQUIRE(levels.has(layer_z - 100));
        REQUIRE(levels.has(layer_z));
        // The outer wall is the first thing printed at either height. The passes differ in width by 0.05 mm.
        const double distance = (levels.first_start.at(layer_z) - levels.first_start.at(layer_z - 100)).norm();
        REQUIRE(distance < 0.1);
    }
}

TEST_CASE("The filament of the split outer wall is printed first on a layer", "[OuterWallPasses]")
{
    // A 30 mm cap on a 10 mm stem, the support of the cap printed with filament 1 and its interface with filament 2.
    // The tool order would start every other layer with the filament the layer below ended with, that is with
    // the interface, and the low passes of the stem would then be printed below it.
    TriangleMesh mesh = make_cube(10, 10, 10);
    mesh.translate(10, 10, 0);
    TriangleMesh cap = make_cube(30, 30, 2);
    cap.translate(0, 0, 10);
    mesh.merge(cap);
    const DynamicPrintConfig config = multifilament_config(2, {
        { "split_outer_wall",                 true },
        { "split_outer_wall_count",           2 },
        { "split_outer_wall_min_width",       "0.4" },
        { "split_outer_wall_max_width",       "0.6" },
        { "split_outer_wall_fallback",        "regular" },
        { "wall_generator",                   "classic" },
        { "layer_height",                     layer_height },
        { "initial_layer_print_height",       layer_height },
        { "elefant_foot_compensation",        0. },
        { "enable_prime_tower",               false },
        { "enable_support",                   true },
        { "support_type",                     "normal(auto)" },
        { "independent_support_layer_height", false },
        { "support_filament",                 1 },
        { "support_interface_filament",       2 },
        { "skirt_loops",                      0 },
        { "brim_type",                        "no_brim" },
    });
    const std::string gcode = Slic3r::Test::slice({ mesh }, config);
    // The interface is printed with the second filament on several layers: the second of them is the first
    // one the tool order would start with the interface.
    size_t interface_changes = 0;
    for (size_t pos = gcode.find("\nT1"); pos != std::string::npos; pos = gcode.find("\nT1", pos + 1))
        ++ interface_changes;
    REQUIRE(interface_changes >= 2);
    const ExtrusionLevels levels = extrusion_levels(gcode);
    REQUIRE(levels.has(5100));
    REQUIRE_FALSE(levels.went_below);
}

TEST_CASE("Walls of holes are split only if enabled", "[OuterWallPasses]")
{
    // Lower passes and top passes of a layer of a cube with a hole.
    auto pass_loops = [](bool holes) {
        DynamicPrintConfig config = passes_config("split_layer");
        config.set_deserialize_strict({ { "split_outer_wall_holes", holes } });
        Print print;
        init_and_process_print({ TestMesh::cube_with_hole }, print, config);
        const Layer &layer = *print.objects().front()->layers()[10];
        REQUIRE(layer.outer_wall_pass_slices);
        std::pair<size_t, size_t> count(0, 0);
        for (const LayerRegion *region : layer.regions()) {
            REQUIRE(region->outer_wall_lower_passes.size() == 1);
            count.first  += count_pass_loops(region->outer_wall_lower_passes.front());
            count.second += count_pass_loops(region->perimeters);
        }
        return count;
    };
    // The contour and the hole.
    REQUIRE(pass_loops(true) == std::pair<size_t, size_t>(2, 2));
    // The contour only.
    REQUIRE(pass_loops(false) == std::pair<size_t, size_t>(1, 1));
}

TEST_CASE("Nothing is split with the option off", "[OuterWallPasses]")
{
    DynamicPrintConfig config = passes_config("split_layer");
    config.set_deserialize_strict({ { "split_outer_wall", false } });
    Print print;
    init_and_process_print({ cube(20.) }, print, config);
    const PrintObject &object = *print.objects().front();
    REQUIRE(object.layers().size() == 100);
    for (const Layer *layer : object.layers()) {
        REQUIRE_FALSE(layer->outer_wall_pass_slices);
        for (const PassPaths &pass : passes_of(*layer))
            REQUIRE(pass.widths.empty());
    }
    const ExtrusionLevels levels = extrusion_levels(print);
    REQUIRE_FALSE(levels.has(300));
    REQUIRE(levels.has(400));
}

TEST_CASE("Passes keep the Z offset of the printer", "[OuterWallPasses]")
{
    DynamicPrintConfig config = passes_config("split_layer");
    config.set_deserialize_strict({ { "z_offset", 0.05 } });
    Print print;
    init_and_process_print({ cube(20.) }, print, config);
    const ExtrusionLevels levels = extrusion_levels(print);
    REQUIRE_FALSE(levels.went_below);
    // The first layer, then the low pass and the rest of the second layer, all of them 0.05 mm higher.
    REQUIRE(levels.has(250));
    REQUIRE(levels.has(350));
    REQUIRE(levels.has(450));
    REQUIRE_FALSE(levels.has(300));
    REQUIRE_FALSE(levels.has(400));
}

TEST_CASE("The top pass is printed before the infill even with the infill first", "[OuterWallPasses]")
{
    DynamicPrintConfig config = passes_config("split_layer");
    config.set_deserialize_strict({ { "is_infill_first", true } });
    Print print;
    init_and_process_print({ cube(20.) }, print, config);
    const ExtrusionLevels levels = extrusion_levels(print);
    for (int layer_z : { 2000, 10400 }) {
        CAPTURE(layer_z);
        REQUIRE(levels.has(layer_z - 100));
        REQUIRE(levels.has(layer_z));
        // The first thing printed at the layer height starts where the low pass did: it is the top pass, not the infill.
        const double distance = (levels.first_start.at(layer_z) - levels.first_start.at(layer_z - 100)).norm();
        REQUIRE(distance < 0.1);
    }
}

TEST_CASE("Nothing is purged into a layer printed in passes", "[OuterWallPasses]")
{
    // The cap on a stem again, its support interface printed with the second filament. The change back to the first
    // filament is purged at a prime tower and, as far as allowed, into the infill, which is set to go first.
    TriangleMesh mesh = make_cube(10, 10, 10);
    mesh.translate(10, 10, 0);
    TriangleMesh cap = make_cube(30, 30, 2);
    cap.translate(0, 0, 10);
    mesh.merge(cap);
    auto slice_cap = [&mesh](bool split) {
        return Slic3r::Test::slice({ mesh }, multifilament_config(2, {
            { "split_outer_wall",                 split },
            { "split_outer_wall_count",           2 },
            { "split_outer_wall_min_width",       "0.4" },
            { "split_outer_wall_max_width",       "0.6" },
            { "split_outer_wall_fallback",        "regular" },
            { "wall_generator",                   "classic" },
            { "layer_height",                     layer_height },
            { "initial_layer_print_height",       layer_height },
            { "elefant_foot_compensation",        0. },
            { "sparse_infill_density",            "15%" },
            { "is_infill_first",                  true },
            { "enable_prime_tower",               true },
            // On the 200 x 200 mm bed of the tests, which the default position is not, and clear of the model.
            { "wipe_tower_x",                     "50" },
            { "wipe_tower_y",                     "50" },
            // One nozzle purged at the tower: what is purged into the object is taken off the tower.
            { "single_extruder_multi_material",   true },
            { "purge_in_prime_tower",             true },
            { "flush_into_infill",                true },
            { "enable_support",                   true },
            { "support_type",                     "normal(auto)" },
            { "independent_support_layer_height", false },
            { "support_filament",                 1 },
            { "support_interface_filament",       2 },
            { "skirt_loops",                      0 },
            { "brim_type",                        "no_brim" },
        }));
    };
    // Without the passes the infill is purged into: the settings do reach the wiping.
    REQUIRE(slice_cap(false).find("; PURGING FINISHED") != std::string::npos);

    const std::string gcode = slice_cap(true);
    // Heights the outer wall was printed at, in microns, and whether an infill was printed at the height of a layer
    // between the low pass and the top pass of its outer wall.
    std::set<int> wall_heights;
    bool          infill_before_top_pass = false;
    // Whether a low pass was the first wall printed after a change to the first filament.
    bool          changed_to_wall_filament = false;
    bool          passes_after_change = false;
    std::string   role;
    GCodeReader   reader;
    reader.parse_buffer(gcode, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        // The role is tagged either way, by the printer.
        for (const std::string tag : { "; FEATURE: ", ";TYPE:" })
            if (line.raw().rfind(tag, 0) == 0)
                role = line.raw().substr(tag.size());
        if (line.raw().rfind("T0", 0) == 0)
            changed_to_wall_filament = true;
        if (! line.extruding(self) || line.dist_XY(self) <= 0.)
            return;
        const int z = int(std::lround(line.new_Z(self) * 1000.));
        if (role.rfind("Outer wall", 0) == 0) {
            // The low pass of a layer lies half a layer below its height.
            if (changed_to_wall_filament && z % 200 == 100)
                passes_after_change = true;
            changed_to_wall_filament = false;
            wall_heights.insert(z);
        } else if (role.find("infill") != std::string::npos && wall_heights.count(z - 100) == 1 && wall_heights.count(z) == 0)
            infill_before_top_pass = true;
    });
    // The stem is printed in passes, a layer started by the change back from the interface included.
    REQUIRE(wall_heights.count(5100) == 1);
    REQUIRE(passes_after_change);
    REQUIRE_FALSE(infill_before_top_pass);
}

TEST_CASE("Changing fuzzy skin plans the layers again", "[OuterWallPasses]")
{
    DynamicPrintConfig config = passes_config("split_layer");
    Print print;
    Model model;
    init_print({ cube(20.) }, print, model, config);
    print.process();
    REQUIRE(print.objects().front()->layers()[1]->outer_wall_pass_slices);

    // Fuzzy skin rules the passes out.
    config.set_deserialize_strict({ { "fuzzy_skin", "external" } });
    print.apply(model, config);
    print.process();
    REQUIRE_FALSE(print.objects().front()->layers()[1]->outer_wall_pass_slices);
    for (const PassPaths &pass : passes_of(*print.objects().front()->layers()[1]))
        REQUIRE(pass.widths.empty());

    config.set_deserialize_strict({ { "fuzzy_skin", "disabled_fuzzy" } });
    print.apply(model, config);
    print.process();
    REQUIRE(print.objects().front()->layers()[1]->outer_wall_pass_slices);
}

TEST_CASE("The support layers are synchronized only where the passes apply", "[OuterWallPasses]")
{
    auto independent_support_layers = [](const char *wall_generator) {
        DynamicPrintConfig config = passes_config("split_layer");
        config.set_deserialize_strict({
            { "wall_generator",                   wall_generator },
            { "enable_support",                   true },
            { "independent_support_layer_height", true },
        });
        Print print;
        Model model;
        init_print({ cube(20.) }, print, model, config);
        return print.config().independent_support_layer_height.value;
    };
    REQUIRE_FALSE(independent_support_layers("classic"));
    // The passes are not generated with Arachne, so the setting is left alone.
    REQUIRE(independent_support_layers("arachne"));
}

TEST_CASE("Settings the passes cannot be printed safely with are rejected", "[OuterWallPasses]")
{
    // What validate() says about one or two cubes with the given changes of the settings.
    auto validate = [](std::initializer_list<ConfigBase::SetDeserializeItem> items, double second_cube_height = 0.) {
        DynamicPrintConfig config = passes_config("split_layer");
        // Relative extrusion needs the extruder reset at every layer to pass the validation at all.
        config.set_deserialize_strict({ { "layer_change_gcode", "G92 E0\n" } });
        config.set_deserialize_strict(items);
        Print print;
        Model model;
        if (second_cube_height > 0.)
            init_print({ cube(20.), make_cube(20., 20., second_cube_height) }, print, model, config);
        else
            init_print({ cube(20.) }, print, model, config);
        return print.validate().string;
    };

    // The same about two cubes with a setting of the second one changed.
    auto validate_second_cube = [](const char *fallback, const ConfigBase::SetDeserializeItem &second_cube_setting) {
        DynamicPrintConfig config = passes_config(fallback);
        config.set_deserialize_strict({ { "layer_change_gcode", "G92 E0\n" }, { "min_layer_height", "0.04" } });
        const std::vector<std::vector<ConfigBase::SetDeserializeItem>> overrides { {}, { second_cube_setting } };
        Print print;
        Model model;
        init_print(std::vector<TriangleMesh>{ cube(20.), cube(20.) }, print, model, config, &overrides);
        return print.validate().string;
    };

    SECTION("the settings of the other tests are accepted") {
        REQUIRE(validate({}).empty());
        REQUIRE(validate({}, 20.).empty());
    }
    SECTION("a skirt around all the objects higher than one layer stands above the low passes") {
        REQUIRE_THAT(validate({ { "skirt_loops", 2 }, { "skirt_height", 3 } }), ContainsSubstring("skirt"));
        REQUIRE_THAT(validate({ { "skirt_loops", 2 }, { "draft_shield", "enabled" } }), ContainsSubstring("draft shield"));
        // One layer of a skirt is below every pass.
        REQUIRE(validate({ { "skirt_loops", 2 }, { "skirt_height", 1 } }).empty());
    }
    SECTION("precise Z height shifts the last layers of objects of different heights") {
        REQUIRE_THAT(validate({ { "precise_z_height", true } }, 20.05), ContainsSubstring("above the lowest pass"));
        REQUIRE(validate({ { "precise_z_height", true } }, 20.).empty());
    }
    SECTION("thin layers another object may split its layers into count as well") {
        // The second cube is printed at half the layer height. Its layers end at the lowest pass of the first
        // cube or at its layer height, but the thin layers it may be split into end in between.
        REQUIRE_THAT(validate_second_cube("split_layer", { "layer_height", 0.1 }), ContainsSubstring("above the lowest pass"));
        REQUIRE(validate_second_cube("regular", { "layer_height", 0.1 }).empty());
    }
    SECTION("a raft under another object has layers of heights of its own") {
        REQUIRE_THAT(validate_second_cube("regular", { "raft_layers", 3 }), ContainsSubstring("raft under another object"));
    }
    SECTION("a width range that makes no printable line") {
        REQUIRE_THAT(validate({ { "split_outer_wall_min_width", "0.01" }, { "split_outer_wall_max_width", "0.01" } }),
                     ContainsSubstring("too small"));
        REQUIRE_THAT(validate({ { "split_outer_wall_max_width", "3" } }), ContainsSubstring("too large"));
    }
}

TEST_CASE("A wall cut for the bridge of a neighbouring region keeps no pass below it", "[OuterWallPasses]")
{
    // Two pillars 10 mm apart and a deck over the gap between them as a part with walls of its own. The deck takes
    // a strip of the pillars to anchor its bridge on, and the walls of the pillars are cut away there.
    TriangleMesh pillars = make_cube(10., 20., 6.);
    TriangleMesh pillar  = make_cube(10., 20., 6.);
    pillar.translate(20.f, 0.f, 0.f);
    pillars.merge(pillar);
    TriangleMesh deck = make_cube(10., 20., 1.);
    deck.translate(10.f, 0.f, 5.f);
    Model model;
    ModelObject *object = model.add_object();
    object->name = "deck.stl";
    object->add_volume(std::move(pillars));
    object->add_volume(std::move(deck))->config.set("wall_loops", 3);
    object->add_instance();
    object->ensure_on_bed();
    Print print;
    print.apply(model, passes_config("regular"));
    print.set_status_silent();
    print.process();

    // Every pass below the layer height has its top pass, and a wall that is no loop any more has neither.
    size_t layers_with_passes = 0;
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions()) {
            const size_t top_passes = count_pass_loops(region->perimeters);
            CAPTURE(layer->print_z, top_passes);
            for (const ExtrusionEntityCollection &pass : region->outer_wall_lower_passes) {
                REQUIRE(pass.entities.size() == top_passes);
                REQUIRE(count_pass_loops(pass) == top_passes);
            }
            if (top_passes > 0)
                ++ layers_with_passes;
        }
    REQUIRE(layers_with_passes > 10);
}
