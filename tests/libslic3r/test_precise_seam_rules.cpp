#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/GCode/PreciseSeamRules.hpp"
#include "libslic3r/GCode/PreciseSeamRulesApply.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::PreciseSeam;

namespace {
using Layers = std::vector<size_t>;
using Code   = RuleDiagnosticCode;

// Every event of the examples is a layer of the object itself.
Layers object_layers(size_t count)
{
    Layers out;
    for (size_t layer = 1; layer <= count; ++layer)
        out.push_back(layer);
    return out;
}

Layers selected(const LayerSelection &selection)
{
    Layers out;
    for (size_t layer = 1; layer < selection.layers.size(); ++layer)
        if (selection.layers[layer])
            out.push_back(layer);
    return out;
}

// Separate name: with select() below, select({}, {}) would pick this overload.
Layers apply_rules(const LayerRules &rules, size_t count = 400) { return selected(select_layers(rules, count, object_layers(count))); }
Layers select(const std::string &name, const std::string &field = {}, size_t count = 400)
{
    return apply_rules(parse_layer_rules(name, field), count);
}

Layers range(size_t first, size_t last, size_t step = 1)
{
    Layers out;
    for (size_t layer = first; layer <= last; layer += step)
        out.push_back(layer);
    return out;
}
Layers join(Layers a, const Layers &b)
{
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

// Rules given only in the field, as one source.
LayerRules field(const std::string &text) { return parse_layer_rules({}, text); }

size_t count_code(const std::vector<RuleDiagnostic> &diagnostics, Code code)
{
    size_t n = 0;
    for (const RuleDiagnostic &d : diagnostics)
        n += d.code == code;
    return n;
}
} // namespace

TEST_CASE("Ranges, parity and modes compose from left to right", "[PreciseSeam][Rules]")
{
    CHECK(select({}, "129-312") == range(129, 312));
    CHECK(select({}, "excl,1-10") == range(11, 400));
    CHECK(select({}, "odd,excl,101-200") == join(range(1, 99, 2), range(201, 399, 2)));
    CHECK(select({}, "excl,even,incl,50") == join(join(range(1, 49, 2), {50}), range(51, 399, 2)));
    CHECK(select({}, "excl,1-400").empty());
}

TEST_CASE("Formulas take the values within the existing layers for x from zero", "[PreciseSeam][Rules]")
{
    CHECK(select({}, "3x+1") == range(1, 400, 3));
    CHECK(select({}, "-2x+3") == Layers{1, 3});
    CHECK(select({}, "x") == range(1, 400));
    CHECK(select({}, "-x+5") == range(1, 5));
    CHECK(select({}, "2x-1") == range(1, 399, 2));
    // Values below 1 for x < 199 are skipped; x = 199 gives layer 1.
    CHECK(select({}, "+2x-397") == range(1, 399, 2));
}

TEST_CASE("Rules without a valid layer set act on every layer", "[PreciseSeam][Rules]")
{
    const LayerRules reversed = field("312-129");
    CHECK(count_code(reversed.diagnostics, Code::InvalidElement) == 1);
    CHECK(apply_rules(reversed) == range(1, 400));

    const LayerRules too_large = field("1x+5000000");
    CHECK(count_code(too_large.diagnostics, Code::NumberOutOfRange) == 1);
    CHECK(apply_rules(too_large) == range(1, 400));

    CHECK(select({}, {}) == range(1, 400));
    CHECK(select("plain name", {}) == range(1, 400));
}

TEST_CASE("A malformed element is discarded without changing the mode", "[PreciseSeam][Rules]")
{
    const LayerRules rules = field("exclbad,1");
    REQUIRE(rules.diagnostics.size() == 1);
    CHECK(rules.diagnostics.front().code == Code::InvalidElement);
    CHECK(rules.diagnostics.front().element == 1);
    CHECK(rules.diagnostics.front().text == "exclbad");
    CHECK(apply_rules(rules) == Layers{1});
}

TEST_CASE("A valid set outside the existing layers is applied and reported", "[PreciseSeam][Rules]")
{
    const LayerRules     rules     = field("500");
    const LayerSelection selection = select_layers(rules, 400, object_layers(400));
    CHECK(selected(selection).empty());
    REQUIRE(selection.diagnostics.size() == 1);
    CHECK(selection.diagnostics.front().code == Code::EmptySet);
}

TEST_CASE("The first type and priority stay in effect over all sources", "[PreciseSeam][Rules]")
{
    const LayerRules blocked = field("red,pr2,1-50");
    CHECK(blocked.type == ModelVolumeType::PRECISE_SEAM_BLOCKED);
    CHECK(blocked.priority == 2);
    CHECK(apply_rules(blocked) == range(1, 50));

    const LayerRules conflict = field("left,right,10");
    CHECK(conflict.type == ModelVolumeType::PRECISE_SEAM_LEFT);
    CHECK(count_code(conflict.diagnostics, Code::TypeConflict) == 1);
    CHECK(apply_rules(conflict) == Layers{10});

    const LayerRules name_first = parse_layer_rules("-^ps=red=ps^-", "green");
    CHECK(name_first.type == ModelVolumeType::PRECISE_SEAM_BLOCKED);
    CHECK(count_code(name_first.diagnostics, Code::TypeConflict) == 1);
    CHECK(name_first.diagnostics.back().source.field);

    // Repeating the same type or priority value is not a conflict.
    const LayerRules repeated = field("gray,grey,pr1,pr01,pr-3");
    CHECK(repeated.type == ModelVolumeType::PRECISE_SEAM_NEUTRAL);
    CHECK(repeated.priority == 1);
    CHECK(count_code(repeated.diagnostics, Code::TypeConflict) == 0);
    CHECK(count_code(repeated.diagnostics, Code::PriorityConflict) == 1);
    CHECK(field("pr-7").priority == -7);
}

TEST_CASE("Every source starts in the include mode while the set carries over", "[PreciseSeam][Rules]")
{
    CHECK(select("-^ps=excl,1-5=ps^-", "20-30") == range(6, 400));
    CHECK(select("-^ps=1-100=ps^-", "excl,40-50") == join(range(1, 39), range(51, 100)));

    const LayerRules two_blocks = parse_layer_rules("-^ps=excl,1=ps^- part -^ps=2=ps^-", {});
    CHECK(count_code(two_blocks.diagnostics, Code::MultipleBlocks) == 1);
    CHECK(apply_rules(two_blocks) == range(2, 400));

    // A type label is not a layer set: the first set, in the field, still picks the initial set.
    const LayerRules label_then_exclusion = parse_layer_rules("-^ps=red=ps^-", "excl,1");
    CHECK(label_then_exclusion.type == ModelVolumeType::PRECISE_SEAM_BLOCKED);
    CHECK(apply_rules(label_then_exclusion) == range(2, 400));
}

TEST_CASE("Name markers are extracted without case and recovered when broken", "[PreciseSeam][Rules]")
{
    CHECK(select("Bracket -^PS= 10 - 20 =Ps^- v2") == range(10, 20));

    const LayerRules reopened = parse_layer_rules("-^ps=odd -^ps=even=ps^-", {});
    CHECK(count_code(reopened.diagnostics, Code::BrokenMarker) == 1);
    CHECK(apply_rules(reopened) == range(2, 400, 2));

    const LayerRules unclosed = parse_layer_rules("a -^ps=1-5", {});
    CHECK(count_code(unclosed.diagnostics, Code::BrokenMarker) == 1);
    CHECK(apply_rules(unclosed) == range(1, 400));

    const LayerRules stray = parse_layer_rules("a =ps^- -^ps=7=ps^-", {});
    CHECK(count_code(stray.diagnostics, Code::BrokenMarker) == 1);
    CHECK(apply_rules(stray) == Layers{7});
}

TEST_CASE("Whitespace is removed and letters are lowered before parsing", "[PreciseSeam][Rules]")
{
    CHECK(select({}, " EXCL ,\t1 2 ") == join(range(1, 11), range(13, 400)));
    CHECK(select({}, "Incl ODD") == range(1, 399, 2));
    CHECK(select({}, ",,5,,") == Layers{5});
}

TEST_CASE("Elements outside the grammar are invalid as a whole", "[PreciseSeam][Rules]")
{
    const std::string text = GENERATE(as<std::string>{}, "0", "0x+5", "3x+-1", "1+2x", "+5", "-5", "pr0", "pr-0", "pr1.", "pr.1",
                                      "inclred", "inclpr1", "exclexcl1", "red1-5", "excl.", "incl.odd", "red.", "1.5", "z",
                                      "100step3", "100-200step0", "100-200step.3", "top0", "top-5", "top.5", "incloff",
                                      "offjunk", "exclnowarn", "nowarn.", "odd.", "1-", "-1-5", "1-5step", "\xD0\xB0",
                                      "z.5-1", "z12.-30", "z1.2.3-4", "z10-20step3", "ztop5", "inclzinbad", "exclzin1-0",
                                      "zinc1", "zz1", "zin", "z1incl", "zin1step2", "z30-12.5", "z-1", "z+1", "z1e2", "z1-",
                                      "incloverride", "overridered", "overrideoff", "relay0", "relay-1", "relay.1", "relay=1",
                                      "relay1-5", "relay1junk", "inclrelay1", "exclrelay1", "relay1red", "relay-", "relayrandom",
                                      "relay-rand", "relay-random0", "relay-random-1", "relay-random.1", "relay-randomjunk",
                                      "relay-random1-5", "inclrelay-random", "relay-random1red", "inclida", "exclidx",
                                      "magnet.5", "magnet2.", "magnet-2", "magnet2-3", "magnetx", "inclmagnet", "exclmagnet2",
                                      "magnet2red", "magnetz2", "magnet1.2.3", "seed7", "relayseed7", "relay1seed7",
                                      "relay-randomseed", "relay-random1seed", "relay-random1seed4294967296",
                                      "relay-random1seed99999999999999999999", "relay-random1seed7x", "relay-random1seed-1",
                                      "relay-random1seed1.5", "relay-randomseed7seed8", "relay-random1seed7red", "modelx",
                                      "bed2", "modelz10", "inclmodel", "exclbed", "bedred", "model1");
    CAPTURE(text);
    const LayerRules rules = field(text);
    REQUIRE(rules.diagnostics.size() == 1);
    CHECK(rules.diagnostics.front().code == Code::InvalidElement);
    CHECK_FALSE(rules.restricts_layers());
    CHECK_FALSE(rules.type.has_value());
    CHECK_FALSE(rules.priority.has_value());
    CHECK_FALSE(rules.relay.has_value());
    CHECK_FALSE(rules.id.has_value());
    CHECK_FALSE(rules.magnet.has_value());
    CHECK_FALSE(rules.off);
    CHECK_FALSE(rules.nowarn);
}

TEST_CASE("Numbers accept leading zeros up to the limit", "[PreciseSeam][Rules]")
{
    CHECK(select({}, "007") == Layers{7});
    CHECK(select({}, "excl") == range(1, 400)); // A mode alone is not a set: all layers remain.
    CHECK(select({}, "excl,excl5") == join(range(1, 4), range(6, 400)));
    CHECK(field("red").type == ModelVolumeType::PRECISE_SEAM_BLOCKED);
    CHECK(field("pr1").priority == 1);
    CHECK(field("pr-2").priority == -2);
    CHECK(select({}, "1-1000000") == range(1, 400));
    CHECK(count_code(field("1000001").diagnostics, Code::NumberOutOfRange) == 1);
}

TEST_CASE("Step, top, all and the flags follow the additional elements", "[PreciseSeam][Rules]")
{
    CHECK(select({}, "100-200step3") == range(100, 199, 3));
    CHECK(select({}, "top5") == range(396, 400));
    CHECK(select({}, "top500") == range(1, 400));
    CHECK(select({}, "excltop5") == range(1, 395));
    CHECK(select({}, "all") == range(1, 400));
    CHECK(select({}, "excl,all").empty());

    const LayerRules off = field("1-5,off");
    CHECK(off.off);
    CHECK(apply_rules(off).empty());
    const LayerRules nowarn = parse_layer_rules("-^ps=nowarn=ps^-", "nowarn");
    CHECK(nowarn.nowarn);
    CHECK(nowarn.diagnostics.empty());
}

TEST_CASE("Top selects the highest layers of the object, not of the plate", "[PreciseSeam][Rules]")
{
    // ByLayer: object A has events 1 and 3 of three, a taller object B adds event 2.
    const LayerRules rules = field("top1");
    CHECK(selected(select_layers(rules, 3, {1, 3})) == Layers{3});
    // A shorter object whose last layer is event 2 of 5 keeps its own top layer.
    CHECK(selected(select_layers(rules, 5, {1, 2})) == Layers{2});
    // An object without numbered layers has nothing to select.
    CHECK(selected(select_layers(rules, 5, {})).empty());
}

TEST_CASE("Support-only events take numbers that parity counts", "[PreciseSeam][Rules]")
{
    // Model layers 0.2 and 0.4 are events 1 and 3; event 2 is support only and has no perimeters.
    const Layers model = {1, 3};
    const auto   on_model = [&](const std::string &text) {
        const LayerSelection selection = select_layers(field(text), 3, model);
        Layers               out;
        for (size_t layer : model)
            if (selection.contains(layer))
                out.push_back(layer);
        return out;
    };
    CHECK(on_model("odd") == Layers{1, 3});
    CHECK(on_model("even").empty());
}

TEST_CASE("No layers means every set is empty", "[PreciseSeam][Rules]")
{
    CHECK(select({}, "all", 0).empty());
    CHECK(select({}, {}, 0).empty());
}

TEST_CASE("Diagnostics name the source and the element position", "[PreciseSeam][Rules]")
{
    const LayerRules rules = parse_layer_rules("-^ps=1,bad=ps^- -^ps=,junk=ps^-", "2,,oops");
    std::vector<RuleDiagnostic> invalid;
    for (const RuleDiagnostic &d : rules.diagnostics)
        if (d.code == Code::InvalidElement)
            invalid.push_back(d);
    REQUIRE(invalid.size() == 3);
    CHECK(invalid[0].source == RuleSource{false, 1});
    CHECK(invalid[0].element == 2);
    CHECK(invalid[1].source == RuleSource{false, 2});
    CHECK(invalid[1].element == 2);
    CHECK(invalid[2].source.field);
    CHECK(invalid[2].element == 3);
    CHECK(invalid[2].text == "oops");
}

TEST_CASE("Rule types regroup the volumes and priorities order each group", "[PreciseSeam][Rules]")
{
    // The object list holds [A Center, B Left, C Blocked]; A is relabelled green without a priority.
    const LayerRules a = field("green"), none = field({});
    const std::vector<RuledModifier> modifiers = {{ModelVolumeType::PRECISE_SEAM_CENTER, &a},
                                                  {ModelVolumeType::PRECISE_SEAM_LEFT, &none},
                                                  {ModelVolumeType::PRECISE_SEAM_BLOCKED, &none}};
    const EffectiveOrder order = effective_order(modifiers);
    CHECK(order.strong == Layers{1});
    CHECK(order.weak == Layers{0, 2});
    CHECK(order.type_changed == Layers{0});
    CHECK(order.reordered.empty());
}

TEST_CASE("Only the prioritised members of reordered pairs are reported", "[PreciseSeam][Rules]")
{
    const LayerRules none = field({}), pr2 = field("pr2"), pr1 = field("pr1");
    const std::vector<RuledModifier> modifiers = {{ModelVolumeType::PRECISE_SEAM_ENFORCED, &none},
                                                  {ModelVolumeType::PRECISE_SEAM_ENFORCED, &pr2},
                                                  {ModelVolumeType::PRECISE_SEAM_ENFORCED, &pr1}};
    const EffectiveOrder order = effective_order(modifiers);
    CHECK(order.weak == Layers{2, 1, 0});
    CHECK(order.reordered == Layers{1, 2});

    // A priority that keeps the tree order changes no pair and is not reported.
    const LayerRules first = field("pr1");
    const EffectiveOrder kept = effective_order({{ModelVolumeType::PRECISE_SEAM_LEFT, &first},
                                                 {ModelVolumeType::PRECISE_SEAM_LEFT, &none}});
    CHECK(kept.strong == Layers{0, 1});
    CHECK(kept.reordered.empty());
}

TEST_CASE("Negative priorities follow the volumes without priority", "[PreciseSeam][Rules]")
{
    const LayerRules low = field("pr-1"), lower = field("pr-5"), none = field({}), high = field("pr3");
    const EffectiveOrder order = effective_order({{ModelVolumeType::PRECISE_SEAM_CENTER, &lower},
                                                  {ModelVolumeType::PRECISE_SEAM_CENTER, &low},
                                                  {ModelVolumeType::PRECISE_SEAM_CENTER, &none},
                                                  {ModelVolumeType::PRECISE_SEAM_CENTER, &high}});
    CHECK(order.strong == Layers{3, 2, 1, 0});
}

TEST_CASE("A weak first priority never passes a strong volume", "[PreciseSeam][Rules]")
{
    const LayerRules weak_first = field("pr1"), strong_last = field("pr-100");
    const EffectiveOrder order = effective_order({{ModelVolumeType::PRECISE_SEAM_ENFORCED, &weak_first},
                                                  {ModelVolumeType::PRECISE_SEAM_CENTER, &strong_last}});
    CHECK(order.strong == Layers{1});
    CHECK(order.weak == Layers{0});
    CHECK(order.reordered.empty());
}

TEST_CASE("Off volumes take no part and nowarn silences tree mismatches", "[PreciseSeam][Rules]")
{
    const LayerRules off = field("off,pr1"), quiet = field("red,pr1,nowarn"), none = field({});
    const EffectiveOrder order = effective_order({{ModelVolumeType::PRECISE_SEAM_NEUTRAL, &none},
                                                  {ModelVolumeType::PRECISE_SEAM_NEUTRAL, &off},
                                                  {ModelVolumeType::PRECISE_SEAM_NEUTRAL, &quiet}});
    CHECK(order.weak == Layers{2, 0});
    CHECK(order.type_changed.empty());
    CHECK(order.reordered.empty());
    // The silenced volume is still known, for the log.
    CHECK(order.suppressed == Layers{2});
}

TEST_CASE("Unit conversion keeps the rules field of every volume", "[PreciseSeam][Rules]")
{
    Model model;
    ModelObject *object = model.add_object();
    object->add_volume(make_cube(10, 10, 10))->name = "part";
    ModelVolume *helper = object->add_volume(make_cube(2, 2, 2));
    helper->set_type(ModelVolumeType::PRECISE_SEAM_CENTER);
    helper->precise_seam_rules = "off";
    ModelVolume *dormant = object->add_volume(make_cube(2, 2, 2));
    dormant->precise_seam_rules = "excl,1-5";
    dormant->set_type(ModelVolumeType::PARAMETER_MODIFIER);
    object->add_instance();
    // Converting only the part recreates every volume of the object.
    ModelObjectPtrs converted;
    object->convert_units(converted, ConversionType::CONV_FROM_INCH, {0});
    REQUIRE(converted.size() == 1);
    REQUIRE(converted.front()->volumes.size() == 3);
    CHECK(converted.front()->volumes[1]->precise_seam_rules == "off");
    CHECK(converted.front()->volumes[2]->precise_seam_rules == "excl,1-5");
    // The model owns and deletes them; ModelObject has no public destructor.
    model.objects.insert(model.objects.end(), converted.begin(), converted.end());
}

TEST_CASE("Only rules text on a Precise Seam volume selects the rules branch", "[PreciseSeam][Rules]")
{
    Model model;
    ModelObject *object = model.add_object();
    object->add_volume(make_cube(10, 10, 10))->name = "part";
    ModelVolume *helper = object->add_volume(make_cube(2, 2, 2));
    helper->set_type(ModelVolumeType::PRECISE_SEAM_CENTER);
    helper->name = "helper odd relay1"; // Rule words outside a block are plain text.
    CHECK_FALSE(uses_precise_seam_rules(*object));
    // A broken block in any letter case still takes the branch, so that it is reported.
    for (const std::string name : {"helper -^PS=odd", "helper odd=Ps^-", "-^ps==ps^-"}) {
        CAPTURE(name);
        helper->name = name;
        CHECK(uses_precise_seam_rules(*object));
    }
    helper->name = "helper";
    helper->precise_seam_rules = " "; // Any field text counts, even one that parses to nothing.
    CHECK(uses_precise_seam_rules(*object));
    // A volume that is not Precise Seam keeps its field dormant.
    helper->set_type(ModelVolumeType::PARAMETER_MODIFIER);
    helper->name = "helper -^ps=odd=ps^-";
    CHECK_FALSE(uses_precise_seam_rules(*object));
}

namespace {
// Model layers 1..count of equal height from `base`, each its own event: band i is (base + (i-1)h, base + ih].
std::vector<LayerBand> grid(size_t count, double height, double base = 0.)
{
    std::vector<LayerBand> out;
    for (size_t i = 1; i <= count; ++i)
        out.push_back({i, base + double(i - 1) * height, base + double(i) * height});
    return out;
}

// The examples of the specification: 400 layers of 0.2 mm, tops 0.2 ... 80.0.
LayerSelection select_heights(const std::string &text, size_t count = 400, double height = 0.2)
{
    return select_layers(field(text), count, object_layers(count), grid(count, height));
}
Layers heights(const std::string &text, size_t count = 400, double height = 0.2) { return selected(select_heights(text, count, height)); }
} // namespace

TEST_CASE("Heights select the layers whose band meets the range", "[PreciseSeam][Rules][Heights]")
{
    CHECK(heights("z0-10") == range(1, 50));
    CHECK(heights("z10-20") == range(50, 100)); // Layer 50, 9.8...10.0, holds the level 10.
    CHECK(heights("z12.5-30") == range(63, 150));
    CHECK(heights("z10-20,excl,odd") == range(50, 100, 2));
    CHECK(heights("excl,z0-1") == range(6, 400));
    CHECK(heights("exclz0-1") == range(6, 400));
    CHECK(heights("z10-20,excl,top5") == range(50, 100));
    CHECK(heights("top5,excl,z79.5-80") == Layers{396, 397});
    // Thinner layers keep the same bounds in mm.
    CHECK(heights("z12.5-30", 800, 0.1) == range(125, 300));
}

TEST_CASE("zin selects only the layers wholly inside the range", "[PreciseSeam][Rules][Heights]")
{
    CHECK(heights("zin10-20") == range(51, 100));
    CHECK(heights("zin12.5-30") == range(64, 150));
    CHECK(heights("inclzin0-1") == range(1, 5));
    // Adjacent ranges share the border layer by intersection only.
    CHECK(heights("z0-10,excl,z10-20") == range(1, 49));
    CHECK(heights("zin0-10,excl,zin10-20") == range(1, 50));
}

TEST_CASE("A single height is the range of one level", "[PreciseSeam][Rules][Heights]")
{
    CHECK(heights("z12.5") == Layers{63});
    CHECK(heights("z12.5-12.5") == Layers{63});
    CHECK(heights("z12.6") == Layers{63}); // A border level belongs to the lower layer.
    const LayerRules single = field("z12.5");
    REQUIRE(single.operations.size() == 1);
    CHECK(single.operations.front().set.first == 12500000);
    CHECK(single.operations.front().set.last == 12500000);
    // No layer fits wholly into one level.
    for (const std::string text : {"zin12.5", "zin12.5-12.5", "zin12.6"}) {
        CAPTURE(text);
        const LayerSelection selection = select_heights(text);
        CHECK(selected(selection).empty());
        CHECK(count_code(selection.diagnostics, Code::EmptySet) == 1);
    }
}

TEST_CASE("Height comparisons lean towards the lower layer by the tolerance", "[PreciseSeam][Rules][Heights]")
{
    CHECK(heights("z10.00005") == Layers{50});
    CHECK(heights("z10.0002") == Layers{51});
    CHECK(heights("z0.00005").empty());
    CHECK(heights("z0.0002") == Layers{1});
    CHECK(heights("zin10.00005-20") == range(51, 100));
    CHECK(heights("zin10.0002-20") == range(52, 100));
    CHECK(heights("zin10-19.99995") == range(51, 100));
    CHECK(heights("zin10-19.9998") == range(51, 99));
    // A band thinner than twice the tolerance passes both inclusion bounds of a range, never of a level.
    const std::vector<LayerBand> thin = {{1, 4.99995, 5.00005}};
    CHECK(selected(select_layers(field("zin5"), 1, {1}, thin)).empty());
    CHECK(selected(select_layers(field("zin4.9999-5.0001"), 1, {1}, thin)) == Layers{1});
}

TEST_CASE("Heights select only numbered model layers of the object", "[PreciseSeam][Rules][Heights]")
{
    // Model layers 0.2 and 0.4 are events 1 and 3; event 2 is support only and has no band.
    const std::vector<LayerBand> model = {{1, 0., 0.2}, {3, 0.2, 0.4}};
    CHECK(selected(select_layers(field("z0-100"), 3, {1, 3}, model)) == Layers{1, 3});
    CHECK(selected(select_layers(field("excl,z0-100"), 3, {1, 3}, model)) == Layers{2});
    // Over a raft and its gap the first model layer starts at 1.0.
    const LayerSelection gap = select_layers(field("z0-0.9"), 10, object_layers(10), grid(10, 0.2, 1.));
    CHECK(selected(gap).empty());
    CHECK(count_code(gap.diagnostics, Code::EmptySet) == 1);
    CHECK(selected(select_layers(field("z0-1.1"), 10, object_layers(10), grid(10, 0.2, 1.))) == Layers{1});
    // Above the object, and without layers at all.
    const LayerSelection above = select_heights("z90-100");
    CHECK(selected(above).empty());
    CHECK(count_code(above.diagnostics, Code::EmptySet) == 1);
    CHECK(selected(select_layers(field("z0-10"), 0, {}, {})).empty());
}

TEST_CASE("Height literals take a point, at most six decimals and the limit", "[PreciseSeam][Rules][Heights]")
{
    CHECK(heights("z012.50-30") == range(63, 150));
    CHECK(heights(" Z 12.5 - 30 ") == range(63, 150));
    CHECK(heights("ZIN10-20") == range(51, 100));
    CHECK(heights("z0.123456-1") == range(1, 5));
    // A decimal comma splits the element: z12 is layer 60, 5-30 are layer numbers.
    CHECK(heights("z12,5-30") == join(range(5, 30), {60}));
    CHECK(field("z1000000").restricts_layers());
    CHECK(field("z0001000000.000000").restricts_layers());
    for (const std::string text : {"z1000000.000001", "z0.1234567-1", "z0.1234560-1", "z1000001", "z0-1000000.5",
                                   "z123456789012345678901234567890"}) {
        CAPTURE(text);
        const LayerRules rules = field(text);
        CHECK(count_code(rules.diagnostics, Code::NumberOutOfRange) == 1);
        CHECK_FALSE(rules.restricts_layers());
    }
}

TEST_CASE("override first in the field replaces the name rules", "[PreciseSeam][Rules][Override]")
{
    const LayerRules replaced = parse_layer_rules("-^ps=red,1-10=ps^-", "override,20-30");
    CHECK(replaced.name_overridden);
    CHECK_FALSE(replaced.type.has_value());
    CHECK(replaced.diagnostics.empty());
    CHECK(apply_rules(replaced) == range(20, 30));

    const LayerRules alone = parse_layer_rules("-^ps=red,1-10=ps^-", "override");
    CHECK(alone.name_overridden);
    CHECK_FALSE(alone.type.has_value());
    CHECK(apply_rules(alone) == range(1, 400));

    // Every effect and every diagnostic of the name goes, broken markers too.
    const LayerRules flags = parse_layer_rules("-^ps=off,pr1,nowarn,green,bad=ps^- -^ps=odd -^ps=even=ps^-", "override");
    CHECK_FALSE(flags.off);
    CHECK_FALSE(flags.nowarn);
    CHECK_FALSE(flags.priority.has_value());
    CHECK_FALSE(flags.type.has_value());
    CHECK(flags.diagnostics.empty());

    // Empty elements and whitespace before it do not count.
    const LayerRules spaced = parse_layer_rules("-^ps=red=ps^-", ", \t , Override ,5");
    CHECK(spaced.name_overridden);
    CHECK_FALSE(spaced.type.has_value());
    CHECK(apply_rules(spaced) == Layers{5});
}

TEST_CASE("The rest of the field after override starts afresh", "[PreciseSeam][Rules][Override]")
{
    CHECK(select("-^ps=1-10=ps^-", "override,excl,5") == join(range(1, 4), range(6, 400)));
    CHECK(parse_layer_rules("-^ps=1-10=ps^-", "override,off").off);
    const LayerSelection empty = select_layers(parse_layer_rules("-^ps=1-10=ps^-", "override,zin12.5"), 400, object_layers(400),
                                               grid(400, 0.2));
    CHECK(selected(empty).empty());
    CHECK(count_code(empty.diagnostics, Code::EmptySet) == 1);
}

TEST_CASE("override elsewhere is an invalid element and keeps the name", "[PreciseSeam][Rules][Override]")
{
    const LayerRules late = parse_layer_rules("-^ps=red=ps^-", "20-30,override");
    CHECK_FALSE(late.name_overridden);
    CHECK(late.type == ModelVolumeType::PRECISE_SEAM_BLOCKED);
    REQUIRE(late.diagnostics.size() == 1);
    CHECK(late.diagnostics.front().code == Code::InvalidElement);
    CHECK(late.diagnostics.front().source.field);
    CHECK(late.diagnostics.front().element == 2);
    CHECK(apply_rules(late) == range(20, 30));

    const LayerRules after_bad = parse_layer_rules("-^ps=red=ps^-", "bad,override,5");
    CHECK_FALSE(after_bad.name_overridden);
    CHECK(after_bad.type == ModelVolumeType::PRECISE_SEAM_BLOCKED);
    CHECK(count_code(after_bad.diagnostics, Code::InvalidElement) == 2);
    CHECK_FALSE(parse_layer_rules("-^ps=red=ps^-", "incl,override").name_overridden);

    const LayerRules in_name = parse_layer_rules("-^ps=override,5=ps^-", {});
    CHECK_FALSE(in_name.name_overridden);
    REQUIRE(in_name.diagnostics.size() == 1);
    CHECK(in_name.diagnostics.front().source == RuleSource{false, 1});
    CHECK(apply_rules(in_name) == Layers{5});
    // The override in the name is not reported once the field replaces the name.
    CHECK(parse_layer_rules("-^ps=override,5=ps^-", "override").diagnostics.empty());

    // A repeat is invalid but leaves the name replaced; element positions count the empty elements.
    const LayerRules repeated = parse_layer_rules("-^ps=red=ps^-", ",override,bad,override");
    CHECK(repeated.name_overridden);
    CHECK_FALSE(repeated.type.has_value());
    REQUIRE(repeated.diagnostics.size() == 2);
    CHECK(repeated.diagnostics[0].element == 3);
    CHECK(repeated.diagnostics[0].text == "bad");
    CHECK(repeated.diagnostics[1].element == 4);
}

TEST_CASE("relay labels a group, and the first label stays in effect", "[PreciseSeam][Rules][Relay]")
{
    CHECK(field("relay").relay == RelayLabel{}); // A group of its own, apart from the numbered ones.
    CHECK(field("relay1").relay == RelayLabel{false, 1});
    CHECK(field("Relay007").relay == RelayLabel{false, 7});
    CHECK(field("relay1000000").relay == RelayLabel{false, 1000000});
    CHECK(count_code(field("relay1000001").diagnostics, Code::NumberOutOfRange) == 1);
    CHECK(field("relay1,relay01").diagnostics.empty());
    for (const std::string text : {"relay1,relay2", "relay1,relay", "relay,relay1"}) {
        CAPTURE(text);
        const LayerRules rules = field(text);
        CHECK(rules.relay == field(text.substr(0, text.find(','))).relay);
        REQUIRE(rules.diagnostics.size() == 1);
        CHECK(rules.diagnostics.front().code == Code::RelayConflict);
        CHECK(rules.diagnostics.front().element == 2);
    }
    // The label carries over from the name to the field, like the type.
    const LayerRules both = parse_layer_rules("-^ps=relay2=ps^-", "relay3,odd");
    CHECK(both.relay == RelayLabel{false, 2});
    CHECK(count_code(both.diagnostics, Code::RelayConflict) == 1);
    CHECK(apply_rules(both) == range(1, 399, 2));
}

TEST_CASE("relay-random labels a group apart from the relay groups", "[PreciseSeam][Rules][Relay]")
{
    CHECK(field("relay-random").relay == RelayLabel{true, 0});
    CHECK(field("Relay - Random 07").relay == RelayLabel{true, 7});
    CHECK(count_code(field("relay-random1000001").diagnostics, Code::NumberOutOfRange) == 1);
    CHECK(field("relay-random1,relay-random01").diagnostics.empty());
    // relay1 and relay-random1 are different groups, so one volume cannot be in both.
    for (const std::string text : {"relay1,relay-random1", "relay-random1,relay1", "relay-random,relay-random1", "relay-random,relay"}) {
        CAPTURE(text);
        const LayerRules rules = field(text);
        CHECK(rules.relay == field(text.substr(0, text.find(','))).relay);
        REQUIRE(rules.diagnostics.size() == 1);
        CHECK(rules.diagnostics.front().code == Code::RelayConflict);
    }
    CHECK(relay_label_text({}) == "relay");
    CHECK(relay_label_text({false, 3}) == "relay3");
    CHECK(relay_label_text({true, 0}) == "relay-random");
    CHECK(relay_label_text({true, 12}) == "relay-random12");
}

TEST_CASE("magnet takes an optional distance in mm, and the first one stays in effect", "[PreciseSeam][Rules][Magnet]")
{
    CHECK(field("magnet").magnet == MagnetRule{});
    CHECK(field("magnet2.5").magnet == MagnetRule{2500000});
    CHECK(field("Magnet 2.5").magnet == MagnetRule{2500000});
    CHECK(field("magnet0.000001").magnet == MagnetRule{1});
    // Zero is a limit of its own, a touch, not the absence of one.
    CHECK(field("magnet0").magnet == MagnetRule{0});
    CHECK(field("magnet0.0").magnet == MagnetRule{0});
    CHECK(field("magnet1000000").magnet == MagnetRule{1000000 * nanometres_per_mm});
    CHECK(count_code(field("magnet1000001").diagnostics, Code::NumberOutOfRange) == 1);
    CHECK(count_code(field("magnet1.0000001").diagnostics, Code::NumberOutOfRange) == 1);
    // Repeats compare as numbers.
    CHECK(field("magnet2.5,magnet02.500000").diagnostics.empty());
    CHECK(field("magnet0,magnet00.000000").diagnostics.empty());
    CHECK(field("magnet,magnet").diagnostics.empty());
    for (const std::string text : {"magnet2,magnet3", "magnet,magnet2", "magnet2,magnet", "magnet0,magnet"}) {
        CAPTURE(text);
        const LayerRules rules = field(text);
        CHECK(rules.magnet == field(text.substr(0, text.find(','))).magnet);
        REQUIRE(rules.diagnostics.size() == 1);
        CHECK(rules.diagnostics.front().code == Code::MagnetConflict);
        CHECK(rules.diagnostics.front().element == 2);
    }
    // An invalid element does not take the place of the first magnet.
    const LayerRules late = field("magnet.5,magnet3");
    CHECK(late.magnet == MagnetRule{3000000});
    CHECK(count_code(late.diagnostics, Code::InvalidElement) == 1);
    CHECK(count_code(late.diagnostics, Code::MagnetConflict) == 0);
    // It keeps the mode and the sets, and carries over from the name to the field; nowarn keeps the conflict.
    CHECK(select({}, "excl,magnet,odd") == range(2, 400, 2));
    const LayerRules both = parse_layer_rules("-^ps=magnet2=ps^-", "magnet3,nowarn");
    CHECK(both.magnet == MagnetRule{2000000});
    REQUIRE(both.diagnostics.size() == 1);
    CHECK(both.diagnostics.front().code == Code::MagnetConflict);
    CHECK(both.diagnostics.front().source == RuleSource{true, 0});
    const LayerRules overridden = parse_layer_rules("-^ps=magnet2=ps^-", "override,magnet3");
    CHECK(overridden.magnet == MagnetRule{3000000});
    CHECK(overridden.diagnostics.empty());
}

TEST_CASE("id names a volume for scripts and changes nothing else", "[PreciseSeam][Rules]")
{
    // Any text may follow id, non-ASCII too, and keeps its case.
    for (const std::string text : {"idFront", "id", "id\xD0\x9F", "idA=ps^-B", "idoff", "idoverride", "idrelay-random1", "id5-1x"}) {
        CAPTURE(text);
        const LayerRules rules = field(text);
        CHECK(rules.diagnostics.empty());
        CHECK(rules.id == text);
        CHECK_FALSE(rules.restricts_layers());
        CHECK_FALSE(rules.type.has_value());
        CHECK_FALSE(rules.relay.has_value());
        CHECK_FALSE(rules.off);
        CHECK_FALSE(rules.name_overridden);
    }
    // The original is trimmed and keeps its inner spaces.
    CHECK(field(" I D Front ,5").id == "I D Front");
    // Repeats compare after normalization: ASCII case and spaces do not count, non-ASCII case does.
    const LayerRules same = field("I D A,ida");
    CHECK(same.diagnostics.empty());
    CHECK(same.id == "I D A");
    for (const std::string text : {"idA,idB", "id\xD0\x9F,id\xD0\xBF", "id,idA"}) {
        CAPTURE(text);
        const LayerRules rules = field(text);
        CHECK(rules.id == text.substr(0, text.find(',')));
        REQUIRE(rules.diagnostics.size() == 1);
        CHECK(rules.diagnostics.front().code == Code::IdConflict);
        CHECK(rules.diagnostics.front().element == 2);
        CHECK(rules.diagnostics.front().text == text.substr(text.find(',') + 1));
    }
    // An id keeps the mode and picks no initial set.
    CHECK(select({}, "excl,idA,odd") == range(2, 400, 2));
    CHECK(select({}, "idA,5") == Layers{5});
    // The id carries over from the name to the field, unless override leaves the name out.
    const LayerRules both = parse_layer_rules("-^ps=idN=ps^-", "idF");
    CHECK(both.id == "idN");
    REQUIRE(both.diagnostics.size() == 1);
    CHECK(both.diagnostics.front().source == RuleSource{true, 0});
    const LayerRules overridden = parse_layer_rules("-^ps=idN=ps^-", "override,idF");
    CHECK(overridden.name_overridden);
    CHECK(overridden.id == "idF");
    CHECK(overridden.diagnostics.empty());
    // override after an id is no longer the first element.
    const LayerRules late = parse_layer_rules("-^ps=5=ps^-", "idF,override");
    CHECK_FALSE(late.name_overridden);
    REQUIRE(late.diagnostics.size() == 1);
    CHECK(late.diagnostics.front().code == Code::InvalidElement);
    CHECK(late.diagnostics.front().element == 2);
    CHECK(apply_rules(late) == Layers{5});
}

TEST_CASE("Name markers and commas bound an id", "[PreciseSeam][Rules]")
{
    // The name block ends at the first closing marker: off falls outside, and the second closing is stray.
    const LayerRules closed = parse_layer_rules("-^ps=idA=ps^-B,off=ps^-", {});
    CHECK(closed.id == "idA");
    CHECK_FALSE(closed.off);
    CHECK(count_code(closed.diagnostics, Code::BrokenMarker) == 1);
    // A second opening discards the block with the id.
    const LayerRules reopened = parse_layer_rules("-^ps=idA-^ps=off=ps^-", {});
    CHECK_FALSE(reopened.id.has_value());
    CHECK(reopened.off);
    CHECK(count_code(reopened.diagnostics, Code::BrokenMarker) == 1);
    // A comma always ends the id.
    const LayerRules split = field("idA,B");
    CHECK(split.id == "idA");
    CHECK(count_code(split.diagnostics, Code::InvalidElement) == 1);
}

namespace {
// One row per layer, one character per member: '1' where the member wins; an empty row has no winner.
std::vector<std::vector<uint8_t>> win_table(const std::vector<std::string> &rows)
{
    std::vector<std::vector<uint8_t>> wins;
    for (const std::string &row : rows) {
        wins.emplace_back();
        for (char c : row)
            wins.back().push_back(c == '1');
    }
    return wins;
}
std::vector<size_t> holders(const std::vector<std::string> &rows, size_t members) { return relay_holders(win_table(rows), members); }
std::vector<size_t> random_holders(const std::vector<std::string> &rows, size_t members, int64_t group)
{
    return relay_random_holders(win_table(rows), members, relay_random_group_seed(group, {}));
}
constexpr size_t nobody = RelayOrder::no_holder;
} // namespace

TEST_CASE("The relay turn passes to the member after each holder", "[PreciseSeam][Rules][Relay]")
{
    CHECK(holders({"111", "111", "111", "111"}, 3) == Layers{0, 1, 2, 0});
    // B cannot on layer 2: A holds again, and the turn stays with B.
    CHECK(holders({"11", "10", "11", "11"}, 2) == Layers{0, 0, 1, 0});
    // Nobody wins on layer 2, or the layer has nothing to win: the turn is kept.
    CHECK(holders({"11", "00", "11", "11"}, 2) == Layers{0, nobody, 1, 0});
    CHECK(holders({"11", "", "11"}, 2) == Layers{0, nobody, 1});
    // The first member cannot first; several in a row cannot.
    CHECK(holders({"01", "11", "11"}, 2) == Layers{1, 0, 1});
    CHECK(holders({"001", "111", "111"}, 3) == Layers{2, 0, 1});
    // B acts only from layer 51: A holds layers 1..50, then they alternate.
    std::vector<std::string> rows(50, "10");
    rows.insert(rows.end(), 4, "11");
    std::vector<size_t> expected(50, 0);
    expected.insert(expected.end(), {1, 0, 1, 0});
    CHECK(holders(rows, 2) == expected);
}

// Expected holders come from an independent model of the hash in SPEC v7 (splitmix64 over the group number and
// the layer index); they pin the hash, so that a change of it or a platform difference shows up here.
TEST_CASE("relay-random gives each layer to a winner other than the previous holder", "[PreciseSeam][Rules][Relay]")
{
    const std::vector<std::string> all(12, "111");
    CHECK(random_holders(all, 3, 0) == Layers{0, 1, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1});
    CHECK(random_holders(all, 3, 1) == Layers{2, 0, 2, 1, 2, 1, 0, 1, 2, 0, 2, 1});
    CHECK(random_holders(all, 3, 2) == Layers{1, 2, 1, 0, 2, 0, 2, 0, 1, 2, 0, 2});
    CHECK(random_holders(all, 3, 1) == random_holders(all, 3, 1));
    // Layers without a winner keep the previous holder, and the hash takes the real layer index. The only
    // winner holds, the previous holder too (layer 8).
    CHECK(random_holders({"111", "", "000", "111", "100", "101", "111", "010", "010", "111"}, 3, 1) ==
          Layers{2, nobody, nobody, 1, 0, 2, 0, 1, 1, 0});
    // Two members alternate after the first choice.
    CHECK(random_holders(std::vector<std::string>(6, "11"), 2, 5) == Layers{1, 0, 1, 0, 1, 0});
    // No holder repeats while another member wins.
    const Layers four = random_holders(std::vector<std::string>(12, "1111"), 4, 3);
    CHECK(four == Layers{3, 0, 1, 3, 2, 1, 0, 2, 1, 0, 1, 2});
    for (size_t layer = 1; layer < four.size(); ++layer)
        CHECK(four[layer] != four[layer - 1]);
}

TEST_CASE("A relay block starts at its holder and wraps around in queue order", "[PreciseSeam][Rules][Relay]")
{
    Model        model;
    ModelObject *object = model.add_object();
    std::vector<const ModelVolume*> v; // A, B, C in the group, then X.
    for (int i = 0; i < 4; ++i)
        v.push_back(object->add_volume(make_cube(2, 2, 2)));
    const ModelVolume *A = v[0], *B = v[1], *C = v[2], *X = v[3];
    RelayOrder order;
    order.entries = {{nullptr, 0}, {X}};
    order.groups  = {{{A, B, C}, {0, 1, 2, nobody}}};
    EffectiveModifiers effective;
    for (const ModelVolume *volume : v)
        effective.emplace(volume, EffectiveModifier{ModelVolumeType::PRECISE_SEAM_CENTER, {}});
    effective.at(B).layers = {true, true, false, true}; // Its rules exclude layer 2.

    using Volumes = std::vector<const ModelVolume*>;
    Volumes out{X}; // Cleared first.
    relay_strong_volumes(order, effective, 0, 2, out);
    CHECK(out == Volumes{A, B, C, X});
    relay_strong_volumes(order, effective, 1, 2, out);
    CHECK(out == Volumes{B, C, A, X});
    // A member outside its layers does not stand in.
    relay_strong_volumes(order, effective, 2, 2, out);
    CHECK(out == Volumes{C, A, X});
    // Without a holder the whole block is left out.
    relay_strong_volumes(order, effective, 3, 2, out);
    CHECK(out == Volumes{X});
    // The entries before the block are its seniors: none here.
    relay_strong_volumes(order, effective, 1, 0, out);
    CHECK(out.empty());
}

TEST_CASE("relay-random takes a seed after the group, outside the group label", "[PreciseSeam][Rules][Relay]")
{
    const LayerRules spaced = field("relay-random 1 seed 7");
    CHECK(spaced.diagnostics.empty());
    CHECK(spaced.relay == RelayLabel{true, 1});
    CHECK(spaced.relay_seed == std::optional<uint32_t>(7));
    CHECK(field("relay-randomseed7").relay == RelayLabel{true, 0});
    CHECK(field("relay-randomseed7").relay_seed == std::optional<uint32_t>(7));
    CHECK(field("relay-random1seed0").relay_seed == std::optional<uint32_t>(0));
    CHECK(field("relay-random1seed4294967295").relay_seed == std::optional<uint32_t>(4294967295u));
    CHECK_FALSE(field("relay-random1").relay_seed.has_value());
    // Repeats compare the values, so leading zeros do not conflict.
    CHECK(field("relay-random01seed007,relay-random1seed7").diagnostics.empty());
    // Another seed, or none, is another element: the first stays in effect with its own seed.
    const LayerRules later = field("relay-random1,relay-random1seed7");
    CHECK(count_code(later.diagnostics, Code::RelayConflict) == 1);
    CHECK_FALSE(later.relay_seed.has_value());
    const LayerRules other = field("relay-random1seed7,relay-random1seed8");
    CHECK(count_code(other.diagnostics, Code::RelayConflict) == 1);
    CHECK(other.relay_seed == std::optional<uint32_t>(7));
    // A malformed seed takes no slot.
    const LayerRules after_error = field("relay-random1seedx,relay-random1seed3");
    CHECK(count_code(after_error.diagnostics, Code::InvalidElement) == 1);
    CHECK(after_error.relay_seed == std::optional<uint32_t>(3));
}

TEST_CASE("The group seed mixes every member seed in queue order", "[PreciseSeam][Rules][Relay]")
{
    using Seeds = std::vector<std::optional<uint32_t>>;
    const std::optional<uint32_t> none;
    // Expected values from _myDocs/_relay/26.10.24_ps-seed-ui/relay_seed_oracle.py.
    CHECK(relay_random_group_seed(0, {}) == 0xE220A8397B1DCDAFull);
    CHECK(relay_random_group_seed(0, Seeds{none, none, none}) == 0xE220A8397B1DCDAFull);
    CHECK(relay_random_group_seed(0, Seeds{7, none, none}) == 0x0F50186E0BA16FA4ull);
    CHECK(relay_random_group_seed(0, Seeds{none, 7, none}) == 0x588ACF71AB60E0A2ull);
    CHECK(relay_random_group_seed(0, Seeds{0, none, none}) == 0xD1C0270687984B37ull);
    CHECK(relay_random_group_seed(0, Seeds{4294967295u, none, none}) == 0x831A81BD93315937ull);
    CHECK(relay_random_group_seed(1, Seeds{none, none, 7}) == 0xCEDF49D26DB06234ull);
    CHECK(relay_random_group_seed(1, Seeds{3, 5, 9}) == 0x379E0E373B6FFD9Full);

    const auto wins = win_table(std::vector<std::string>(12, "111"));
    const auto run  = [&wins](const Seeds &seeds) { return relay_random_holders(wins, 3, relay_random_group_seed(0, seeds)); };
    // Without seeds the holders stay those of the group number alone.
    CHECK(run(Seeds{none, none, none}) == Layers{0, 1, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1});
    CHECK(run(Seeds{7, none, none}) == Layers{1, 0, 1, 0, 2, 0, 2, 0, 1, 2, 0, 2});
    // The same seed on another member, and seed0 against no seed, give other holders here.
    CHECK(run(Seeds{none, 7, none}) == Layers{2, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0, 1});
    CHECK(run(Seeds{0, none, none}) == Layers{0, 2, 1, 0, 2, 1, 0, 1, 2, 1, 0, 1});
    // An only winner fixes the phase of two members, whatever the seed.
    std::vector<std::string> forced(8, "11");
    forced.front() = "10";
    for (uint32_t seed : {0u, 1u, 7u, 4294967295u})
        CHECK(relay_random_holders(win_table(forced), 2, relay_random_group_seed(0, Seeds{seed, none})) == Layers{0, 1, 0, 1, 0, 1, 0, 1});
}

TEST_CASE("model counts the model layers through their numbers", "[PreciseSeam][Rules][Origin]")
{
    // 8 numbers; model layer 2 has none, model layer 3 holds number 6.
    const Layers model_to_bed{3, 0, 6};
    const auto   evaluate = [&](const std::string &text) { return select_layers(field(text), 8, {3, 6}, {}, model_to_bed); };
    const LayerSelection unnumbered = evaluate("model,2");
    CHECK(selected(unnumbered).empty());
    CHECK(count_code(unnumbered.diagnostics, Code::EmptySet) == 1);
    CHECK(selected(evaluate("model,3")) == Layers{6});
    CHECK(selected(evaluate("model,all")) == Layers{3, 6});
    CHECK(selected(evaluate("model,odd")) == Layers{3, 6});
    CHECK(selected(evaluate("model,1,bed,6,excl,model,1")) == Layers{6});
    // The bed universe stays: exclusion starts from every number.
    CHECK(selected(evaluate("excl,model,all")) == Layers{1, 2, 4, 5, 7, 8});
    CHECK(selected(evaluate("bed,4")) == Layers{4});
    // The mode carries over the origin and the origin over the mode.
    CHECK(selected(evaluate("model,excl,odd,bed,2")) == Layers{1, 4, 5, 7, 8});
    // An element all of whose layers were already set is not empty.
    CHECK(evaluate("model,3,bed,6").diagnostics.empty());
    // topK takes the highest numbered model layers in either origin.
    CHECK(selected(evaluate("model,top2")) == Layers{3, 6});
    CHECK(selected(select_layers(field("model,top1"), 8, {3, 6}, {}, Layers{3, 6, 0})) == Layers{6});
    // Without a raft the numbers may still differ: other objects take numbers.
    CHECK(selected(select_layers(field("bed,2"), 5, {1, 3, 5}, {}, Layers{1, 3, 5})) == Layers{2});
    CHECK(selected(select_layers(field("model,2"), 5, {1, 3, 5}, {}, Layers{1, 3, 5})) == Layers{3});
}

TEST_CASE("Every source starts from the bed, and the origin changes only the sets", "[PreciseSeam][Rules][Origin]")
{
    const Layers model_to_bed{4, 5, 6, 7, 8}; // A raft takes numbers 1..3.
    const auto   evaluate = [&](const LayerRules &rules) { return selected(select_layers(rules, 8, {4, 5, 6, 7, 8}, {}, model_to_bed)); };
    CHECK(evaluate(parse_layer_rules("a -^ps=model,1=ps^- -^ps=2=ps^-", "3")) == Layers{2, 3, 4});
    CHECK(evaluate(parse_layer_rules("a -^ps=model=ps^-", "1")) == Layers{1});
    // A malformed element keeps the origin.
    CHECK(evaluate(field("model,modelx,2")) == Layers{5});
    // Only flags restrict nothing.
    CHECK_FALSE(field("model").restricts_layers());
    CHECK(field("model,bed,model").diagnostics.empty());
    // Other elements do not count from the model.
    const LayerRules others = field("model,relay-random1seed7,magnet2.5,pr1,odd");
    CHECK(others.diagnostics.empty());
    CHECK(others.relay_seed == std::optional<uint32_t>(7));
    CHECK(others.magnet == MagnetRule{2500000});
    CHECK(others.priority == std::optional<int>(1));
    CHECK(evaluate(others) == Layers{4, 6, 8});
}

TEST_CASE("model heights start at the bottom of the first model layer", "[PreciseSeam][Rules][Origin][Heights]")
{
    // A raft and its gap end at 0.7; the first model layer is 0.3 mm, the next ones 0.2 mm.
    const std::vector<LayerBand> bands{{4, 0.7, 1.0}, {5, 1.0, 1.2}, {6, 1.2, 1.4}, {7, 1.4, 1.6}, {8, 1.6, 1.8}};
    const Layers model_to_bed{4, 5, 6, 7, 8};
    const auto   evaluate = [&](const std::string &text) { return select_layers(field(text), 8, {4, 5, 6, 7, 8}, bands, model_to_bed, 0.7); };
    CHECK(selected(evaluate("model,z0-0.3")) == Layers{4});
    CHECK(selected(evaluate("model,zin0.3-0.7")) == Layers{5, 6});
    CHECK(selected(evaluate("model,z0.5")) == Layers{5});
    CHECK(selected(evaluate("model,z0.50005")) == Layers{5});
    // The level of the bottom has no layer below it, and a single height holds no layer.
    CHECK(selected(evaluate("model,z0")).empty());
    CHECK(selected(evaluate("model,zin0")).empty());
    // From the bed the same range lies under the raft top.
    const LayerSelection bed = evaluate("z0-0.3");
    CHECK(selected(bed).empty());
    CHECK(count_code(bed.diagnostics, Code::EmptySet) == 1);
    CHECK(selected(evaluate("model,z0-0.3,bed,z1.7")) == Layers{4, 8});
    CHECK(selected(evaluate("excl,model,z0-0.5")) == Layers{1, 2, 3, 6, 7, 8});
    // A model bottom of zero makes both origins measure alike.
    CHECK(selected(select_layers(field("model,z0.1-0.3"), 3, {1, 2, 3}, grid(3, 0.2), {1, 2, 3}, 0.)) ==
          selected(select_layers(field("z0.1-0.3"), 3, {1, 2, 3}, grid(3, 0.2), {1, 2, 3}, 0.)));
}
