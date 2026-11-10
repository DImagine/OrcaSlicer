#include <catch2/catch_all.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "test_helpers.hpp"

using namespace Slic3r;
using namespace Slic3r::Test;
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
