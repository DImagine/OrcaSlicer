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
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <cmath>
#include <cstddef>
#include <set>
#include <string>
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
    const std::string gcode = Slic3r::Test::gcode(print);

    // Heights the nozzle extrudes at, and whether it ever extrudes below something already printed.
    std::set<int> extrusion_z_microns;
    double        highest_z     = 0.;
    bool          went_below    = false;
    GCodeReader   reader;
    reader.apply_config(print.config());
    reader.parse_buffer(gcode, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        if (line.extruding(self) && line.dist_XY(self) > 0.) {
            const double z = line.new_Z(self);
            extrusion_z_microns.insert(int(std::lround(z * 1000.)));
            if (z < highest_z - 1e-3)
                went_below = true;
            if (z > highest_z)
                highest_z = z;
        }
    });

    // The nozzle never returns below a height it already printed at.
    REQUIRE_FALSE(went_below);
    // The first layer is printed at its height only, the layers above also half a layer below.
    REQUIRE(extrusion_z_microns.count(100) == 0);
    REQUIRE(extrusion_z_microns.count(200) == 1);
    REQUIRE(extrusion_z_microns.count(300) == 1);
    REQUIRE(extrusion_z_microns.count(400) == 1);
    REQUIRE(extrusion_z_microns.count(19900) == 1);
    REQUIRE(extrusion_z_microns.count(20000) == 1);
}
