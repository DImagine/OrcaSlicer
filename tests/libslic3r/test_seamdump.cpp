#include <catch2/catch_test_macros.hpp>
#include "libslic3r/GCode/SeamDump.hpp"
#include "libslic3r/GCode/SeamDumpPrint.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"

#include <boost/filesystem.hpp>
#include <boost/nowide/fstream.hpp>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::PreciseSeam;

namespace {

Polygon polygon(const std::vector<Vec2d> &points)
{
    Polygon out;
    for (const Vec2d &p : points)
        out.points.push_back(Point::new_scale(p));
    return out;
}

Polygon rectangle(double x0, double y0, double x1, double y1) { return polygon({{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}); }

Polygon circle(double cx, double cy, double r, size_t n)
{
    std::vector<Vec2d> points;
    for (size_t i = 0; i < n; ++i)
        points.emplace_back(cx + r * std::cos(2. * PI * double(i) / double(n)), cy + r * std::sin(2. * PI * double(i) / double(n)));
    return polygon(points);
}

DumpLayer layer(size_t number, std::vector<DumpLoop> loops)
{
    DumpLayer out;
    out.number  = number;
    out.print_z = 0.2 * double(number);
    out.slice_z = out.print_z - 0.1;
    out.height  = 0.2;
    out.loops   = std::move(loops);
    return out;
}

DumpObject object(std::vector<DumpLayer> layers)
{
    DumpObject out;
    out.name   = "Part";
    out.layers = std::move(layers);
    return out;
}

std::string dump(const std::vector<DumpObject> &objects, const DumpSettings &settings = {})
{
    const DumpResult result = write_seamdump(objects, settings);
    REQUIRE_FALSE(result.svg.empty());
    return result.svg;
}

// The value of the first attribute `name` at or after the element with id `id` (the root for an empty id).
std::string attribute(const std::string &svg, const std::string &id, const std::string &name)
{
    const size_t from = id.empty() ? 0 : svg.find("id=\"" + id + "\"");
    REQUIRE(from != std::string::npos);
    const size_t end = svg.find('>', from);
    const size_t at  = svg.find(" " + name + "=\"", from);
    if (at == std::string::npos || at > end)
        return "<none>";
    const size_t begin = at + name.size() + 3;
    return svg.substr(begin, svg.find('"', begin) - begin);
}

} // namespace

TEST_CASE("seamdump numbers have at most two decimals and no trailing zeros", "[PreciseSeam][Seamdump]")
{
    CHECK(seamdump_number(0) == "0");
    CHECK(seamdump_number(1210) == "12.1");
    CHECK(seamdump_number(4025) == "40.25");
    CHECK(seamdump_number(500) == "5");
    CHECK(seamdump_number(-5) == "-0.05");
    CHECK(seamdump_number(-100) == "-1");
}

TEST_CASE("A loop starts at its smallest point, counterclockwise when outer", "[PreciseSeam][Seamdump]")
{
    // Clockwise from the top right corner, with points along the edges for the simplification to drop.
    const Polygon square = polygon({{10, 10}, {10, 5}, {10, 0}, {5, 0}, {0, 0}, {0, 5}, {0, 10}, {5, 10}});
    const std::string svg = dump({object({layer(1, {{square}})})});
    CHECK(attribute(svg, "O1-Z1-L1", "d") == "M0 0 L10 0 L10 10 L0 10 Z");
    CHECK(attribute(svg, "O1-Z1-L1", "data-kind") == "outer");
    CHECK(attribute(svg, "O1-Z1-L1", "data-source-loop") == "0");
    CHECK(attribute(svg, "O1-Z1-L1", "data-len") == "40");
    CHECK(attribute(svg, "O1-Z1-L1", "data-c") == "0,0,90,1;10,0,90,1;10,10,90,1;0,10,90,1");
    CHECK(attribute(svg, "O1-Z1", "data-z") == "0.2");
    // The same loop given counterclockwise from another point is written the same.
    const std::string again = dump({object({layer(1, {{rectangle(0, 0, 10, 10)}})})});
    CHECK(attribute(again, "O1-Z1-L1", "d") == "M0 0 L10 0 L10 10 L0 10 Z");
    // Fewer candidates keep the sharpest corners nearest to M, in path order.
    DumpSettings two;
    two.candidates = 2;
    CHECK(attribute(dump({object({layer(1, {{square}})})}, two), "O1-Z1-L1", "data-c") == "0,0,90,1;10,0,90,1");
    DumpSettings none;
    none.candidates = 0;
    CHECK(attribute(dump({object({layer(1, {{square}})})}, none), "O1-Z1-L1", "data-c") == "<none>");
}

TEST_CASE("A hole runs clockwise, names its outer loop and has reflex corners", "[PreciseSeam][Seamdump]")
{
    DumpLoop hole{rectangle(5, 5, 15, 15), true};
    const std::string svg = dump({object({layer(1, {{rectangle(0, 0, 20, 20)}, hole})})});
    CHECK(attribute(svg, "O1-Z1-L1", "data-kind") == "outer");
    CHECK(attribute(svg, "O1-Z1-L2", "data-kind") == "hole");
    CHECK(attribute(svg, "O1-Z1-L2", "data-parent") == "L1");
    CHECK(attribute(svg, "O1-Z1-L2", "d") == "M5 5 L5 15 L15 15 L15 5 Z");
    CHECK(attribute(svg, "O1-Z1-L2", "data-c") == "5,5,270,0;5,15,270,0;15,15,270,0;15,5,270,0");
}

TEST_CASE("Round loops are written with arcs within the tolerance", "[PreciseSeam][Seamdump]")
{
    const std::string svg = dump({object({layer(1, {{circle(20, 20, 5, 720)}})})});
    const std::string d   = attribute(svg, "O1-Z1-L1", "d");
    CHECK(d.find(" A") != std::string::npos);
    CHECK(d.find(" C") == std::string::npos);
    CHECK(attribute(svg, "O1-Z1-L1", "data-arcs") != "<none>");
    CHECK(std::abs(std::stod(attribute(svg, "O1-Z1-L1", "data-len")) - 2. * PI * 5.) < 0.1);
    // A smooth loop has no corners.
    CHECK(attribute(svg, "O1-Z1-L1", "data-c") == "<none>");
}

TEST_CASE("Loop numbers carry over one to one and start afresh at a split", "[PreciseSeam][Seamdump]")
{
    const DumpLoop whole{rectangle(0, 0, 20, 10)}, left{rectangle(0, 0, 8, 10)}, right{rectangle(12, 0, 20, 10)};
    const DumpLoop moved{rectangle(1, 0, 21, 10)};
    const std::string svg = dump({object({layer(1, {whole}), layer(2, {moved}), layer(3, {right, left}), layer(4, {}), layer(5, {left})})});
    CHECK(attribute(svg, "O1-Z2-L1", "data-from") == "<none>");
    // The split gives new numbers in the order of their points, with the loop they came from.
    CHECK(attribute(svg, "O1-Z3-L2", "d") == "M0 0 L8 0 L8 10 L0 10 Z");
    CHECK(attribute(svg, "O1-Z3-L2", "data-from") == "L1");
    CHECK(attribute(svg, "O1-Z3-L3", "data-from") == "L1");
    CHECK(attribute(svg, "O1-Z3-L3", "data-source-loop") == "0");
    // No number carries over an empty layer.
    CHECK(attribute(svg, "O1-Z5-L4", "data-from") == "<none>");
    // Holes and outer loops share the numbers of an object.
    const DumpLoop hole{rectangle(5, 2, 15, 8), true};
    const std::string both = dump({object({layer(1, {whole, hole}), layer(2, {whole, hole})})});
    CHECK(both.find("id=\"O1-Z1-L1\"") != std::string::npos);
    CHECK(both.find("id=\"O1-Z1-L2\"") != std::string::npos);
}

TEST_CASE("A layer equal to the last one written refers to it", "[PreciseSeam][Seamdump]")
{
    const DumpLoop square{rectangle(0, 0, 10, 10)};
    const std::string svg = dump({object({layer(1, {square}), layer(2, {square}), layer(3, {square}), layer(4, {{rectangle(0, 0, 10, 11)}})})});
    CHECK(attribute(svg, "O1-Z2", "data-same") == "O1-Z1");
    CHECK(attribute(svg, "O1-Z3", "data-same") == "O1-Z1");
    CHECK(attribute(svg, "O1-Z4", "data-same") == "<none>");
    CHECK(svg.find("id=\"O1-Z2-L1\"") == std::string::npos);
}

TEST_CASE("Paint is written along the published path and marks the candidates", "[PreciseSeam][Seamdump]")
{
    DumpLoop square{rectangle(0, 0, 10, 10)};
    square.paint = {{5., 15., false}, {32., 34., true}};
    const std::string svg = dump({object({layer(1, {square})})});
    CHECK(attribute(svg, "O1-Z1-L1", "data-paint") == "enforce:5-15;block:32-34");
    CHECK(attribute(svg, "O1-Z1-L1", "data-c") == "0,0,90,1;10,0,90,1,e;10,10,90,1;0,10,90,1");
    CHECK(svg.find("class=\"enforce\"") != std::string::npos);
    CHECK(svg.find("class=\"block\"") != std::string::npos);
    // Paint across the start is split there; blocked paint wins over enforced.
    DumpLoop wrapped{rectangle(0, 0, 10, 10)};
    wrapped.paint = {{38., 42., false}, {0., 40., true}};
    CHECK(attribute(dump({object({layer(1, {wrapped})})}), "O1-Z1-L1", "data-paint") == "block:0-40");
}

TEST_CASE("The header counts its own bytes and names the slicing", "[PreciseSeam][Seamdump]")
{
    const DumpLoop square{rectangle(0, 0, 10, 10)};
    const std::string svg = dump({object({layer(1, {square})})});
    CHECK(attribute(svg, "", "data-bytes") == std::to_string(svg.size()));
    CHECK(attribute(svg, "", "data-token-estimate") == std::to_string(svg.size() / 4));
    CHECK(attribute(svg, "", "data-dump-id").size() == 64);
    CHECK(attribute(svg, "", "data-tolerance") == "0.1");
    CHECK(attribute(svg, "", "data-candidates-max") == "8");
    CHECK(attribute(svg, "", "data-layer-origin") == "model");
    // Other layer heights with the same loops: another fingerprint.
    DumpObject higher = object({layer(1, {square})});
    higher.layers.front().print_z += 0.05;
    const std::string other = dump({higher});
    CHECK(attribute(other, "", "data-fingerprint") != attribute(svg, "", "data-fingerprint"));
    // The same input gives the same bytes.
    CHECK(dump({object({layer(1, {square})})}) == svg);
    // Names are escaped.
    DumpObject named = object({layer(1, {square})});
    named.name       = "a<b>&\"c\"";
    CHECK(attribute(dump({named}), "O1", "data-object") == "a&lt;b&gt;&amp;&quot;c&quot;");
}

TEST_CASE("Above the size limit the file is built again coarser, or not at all", "[PreciseSeam][Seamdump]")
{
    const std::vector<DumpObject> objects{object({layer(1, {{circle(20, 20, 5, 720)}}), layer(2, {{circle(20, 20, 6, 720)}})})};
    const std::string             full = dump(objects);
    DumpSettings                  enough;
    enough.limit = full.size();
    const DumpResult fits = write_seamdump(objects, enough);
    CHECK_FALSE(fits.svg.empty());
    CHECK(fits.tolerance == 0.1);
    DumpSettings tight;
    tight.limit = 100;
    const DumpResult none = write_seamdump(objects, tight);
    CHECK(none.svg.empty());
    CHECK(none.size > 100);
    // A cancel stops the work.
    DumpSettings cancelled;
    cancelled.throw_if_canceled = [] { throw std::runtime_error("cancelled"); };
    CHECK_THROWS(write_seamdump(objects, cancelled));
}

TEST_CASE("Loops that touch themselves or each other are left out", "[PreciseSeam][Seamdump]")
{
    // Two squares meeting in one point, as one loop.
    const Polygon eight = polygon({{0, 0}, {2, 0}, {2, 2}, {0, 2}, {0, 0}, {-2, 0}, {-2, -2}, {0, -2}});
    const std::string svg = dump({object({layer(1, {{eight}, {rectangle(10, 0, 12, 2)}})})});
    CHECK(attribute(svg, "", "data-dropped") == "1");
    CHECK(attribute(svg, "", "data-dropped-loops") == "O1-Z1-S0:topology");
    CHECK(attribute(svg, "O1-Z1-L2", "data-source-loop") == "1");
    // Two loops along one side.
    const std::string side = dump({object({layer(1, {{rectangle(0, 0, 5, 5)}, {rectangle(5, 0, 10, 5)}, {rectangle(20, 0, 22, 2)}})})});
    CHECK(attribute(side, "", "data-dropped") == "2");
    // A hole that the rounding would move out of its outer loop: here it touches, and both go.
    const std::string near = dump({object({layer(1, {{rectangle(0, 0, 10, 10)}, {rectangle(0.001, 2, 3, 3), true}})})});
    CHECK(attribute(near, "", "data-dropped") != "0");
}

TEST_CASE("Corners turn at least 5 degrees", "[PreciseSeam][Seamdump]")
{
    // A bottom edge bent at its middle by about 4.6 or 6.9 degrees.
    const std::string slight = dump({object({layer(1, {{polygon({{0, 0}, {5, -0.2}, {10, 0}, {10, 10}, {0, 10}})}})})});
    CHECK(attribute(slight, "O1-Z1-L1", "data-c").find("5,-0.2,") == std::string::npos);
    const std::string clear = dump({object({layer(1, {{polygon({{0, 0}, {5, -0.3}, {10, 0}, {10, 10}, {0, 10}})}})})});
    CHECK(attribute(clear, "O1-Z1-L1", "data-c").find("5,-0.3,") != std::string::npos);
}

TEST_CASE("Paint that ends at M marks the candidate there", "[PreciseSeam][Seamdump]")
{
    DumpLoop square{rectangle(0, 0, 10, 10)};
    square.paint = {{38., 40., false}};
    const std::string svg = dump({object({layer(1, {square})})});
    CHECK(attribute(svg, "O1-Z1-L1", "data-paint") == "enforce:38-40");
    CHECK(attribute(svg, "O1-Z1-L1", "data-c").rfind("0,0,90,1,e;", 0) == 0);
}

TEST_CASE("The fingerprint records the paint source and the reach", "[PreciseSeam][Seamdump]")
{
    const DumpLoop square{rectangle(0, 0, 10, 10)};
    const std::string base = attribute(dump({object({layer(1, {square})})}), "", "data-fingerprint");
    DumpObject painted     = object({layer(1, {square})});
    painted.paint_source   = "facets";
    CHECK(attribute(dump({painted}), "", "data-fingerprint") != base);
    DumpObject reach = object({layer(1, {square})});
    reach.layers.front().loops.front().reach = 0.2;
    CHECK(attribute(dump({reach}), "", "data-fingerprint") != base);
}

TEST_CASE("Object file names are portable and numbered when equal", "[PreciseSeam][Seamdump]")
{
    CHECK(seamdump_object_names({"Boat"}) == std::vector<std::string>{"Boat"});
    // Equal names, also without letter case, are all numbered in list order.
    CHECK(seamdump_object_names({"Part", "Box", "part", "Part"}) == std::vector<std::string>{"Part (1)", "Box", "part (2)", "Part (3)"});
    // Forbidden characters, trailing dots and spaces, reserved and empty names.
    CHECK(seamdump_object_names({"A/B", "A\\B"}) == std::vector<std::string>{"A_B (1)", "A_B (2)"});
    CHECK(seamdump_object_names({"x. ", "CON", "", "a:b?"}) == std::vector<std::string>{"x", "_CON", "_", "a_b_"});
    // A name already taken by another object gets the next free number.
    CHECK(seamdump_object_names({"Name", "Name (1)", "Name"}) == std::vector<std::string>{"Name (1)", "Name (1) (1)", "Name (2)"});
    // Letter case beyond ASCII counts as file systems without case see it.
    CHECK(seamdump_fold("\xD0\x94\xD0\xB5\xD1\x82\xD0\xB0\xD0\xBB\xD1\x8C") == "\xD0\xB4\xD0\xB5\xD1\x82\xD0\xB0\xD0\xBB\xD1\x8C");
    CHECK(seamdump_object_names({"\xD0\x94\xD0\xB5\xD1\x82\xD0\xB0\xD0\xBB\xD1\x8C", "\xD0\xB4\xD0\xB5\xD1\x82\xD0\xB0\xD0\xBB\xD1\x8C"}) ==
          std::vector<std::string>{"\xD0\x94\xD0\xB5\xD1\x82\xD0\xB0\xD0\xBB\xD1\x8C (1)", "\xD0\xB4\xD0\xB5\xD1\x82\xD0\xB0\xD0\xBB\xD1\x8C (2)"});
}

TEST_CASE("seamdump file names follow the request and the project", "[PreciseSeam][Seamdump]")
{
    SeamdumpJob plate;
    plate.files.push_back({{}, "<svg/>"});
    CHECK(seamdump_file_names(plate, "Project") == std::vector<std::string>{"Project.seamdump.svg"});
    plate.request.file = "Boat, v2.svg";
    CHECK(seamdump_file_names(plate, "Project") == std::vector<std::string>{"Boat, v2.svg"});
    SeamdumpJob each;
    each.request.each = true;
    each.request.file = "boat.SVG";
    each.files        = {{"Hull (1)", "1"}, {"Cabin", "2"}};
    CHECK(seamdump_file_names(each, "Project") == std::vector<std::string>{"boat-Hull (1).svg", "boat-Cabin.svg"});
    each.request.file.reset();
    CHECK(seamdump_file_names(each, "Project") == std::vector<std::string>{"Project.seamdump-Hull (1).svg", "Project.seamdump-Cabin.svg"});
}

TEST_CASE("A seamdump file replaces the old one only when complete", "[PreciseSeam][Seamdump]")
{
    const boost::filesystem::path folder = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("seamdump-%%%%-%%%%");
    boost::filesystem::create_directories(folder);
    const auto read = [&folder](const std::string &name) {
        boost::nowide::ifstream in((folder / name).string(), std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };
    const boost::filesystem::path target = folder / "Boat, v2.svg";
    CHECK(write_seamdump_file(target, "<svg/>"));
    CHECK(read("Boat, v2.svg") == "<svg/>");
    // The same file again: replaced.
    CHECK(write_seamdump_file(target, "<svg>2</svg>"));
    CHECK(read("Boat, v2.svg") == "<svg>2</svg>");
    // A write that cannot happen leaves the old file and no temporary one.
    CHECK_FALSE(write_seamdump_file(folder / "missing" / "Boat, v2.svg", "<svg>3</svg>"));
    CHECK(read("Boat, v2.svg") == "<svg>2</svg>");
    for (const auto &entry : boost::filesystem::directory_iterator(folder))
        CHECK(entry.path().extension() != ".tmp");
    boost::filesystem::remove_all(folder);
}
