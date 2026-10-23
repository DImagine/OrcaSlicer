#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/GCode/PreciseSeamRules.hpp"
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
    CHECK(select({}, "Incl. ODD.") == range(1, 399, 2));
    CHECK(select({}, ",,5,,") == Layers{5});
}

TEST_CASE("Elements outside the grammar are invalid as a whole", "[PreciseSeam][Rules]")
{
    const std::string text = GENERATE(as<std::string>{}, "0", "0x+5", "3x+-1", "1+2x", "+5", "-5", "pr0", "pr-0", "pr1.", "pr..1",
                                      "inclred", "inclpr1", "exclexcl1", "red1-5", "incl..", "1.5", "z10", "100step3",
                                      "100-200step0", "top0", "top-5", "top5.", "incloff", "offjunk", "exclnowarn", "odd..",
                                      "1-", "-1-5", "1-5step", "\xD0\xB0");
    CAPTURE(text);
    const LayerRules rules = field(text);
    REQUIRE(rules.diagnostics.size() == 1);
    CHECK(rules.diagnostics.front().code == Code::InvalidElement);
    CHECK_FALSE(rules.restricts_layers());
    CHECK_FALSE(rules.type.has_value());
    CHECK_FALSE(rules.priority.has_value());
    CHECK_FALSE(rules.off);
    CHECK_FALSE(rules.nowarn);
}

TEST_CASE("Keywords accept one dot and numbers accept leading zeros", "[PreciseSeam][Rules]")
{
    CHECK(select({}, "007") == Layers{7});
    CHECK(select({}, "excl.") == range(1, 400)); // A mode alone is not a set: all layers remain.
    CHECK(select({}, "excl.,excl.5") == join(range(1, 4), range(6, 400)));
    CHECK(field("red.").type == ModelVolumeType::PRECISE_SEAM_BLOCKED);
    CHECK(field("pr.1").priority == 1);
    CHECK(field("pr.-2").priority == -2);
    CHECK(select({}, "1-1000000") == range(1, 400));
    CHECK(count_code(field("1000001").diagnostics, Code::NumberOutOfRange) == 1);
}

TEST_CASE("Step, top, all and the flags follow the additional elements", "[PreciseSeam][Rules]")
{
    CHECK(select({}, "100-200step3") == range(100, 199, 3));
    CHECK(select({}, "100-200step.3") == range(100, 199, 3));
    CHECK(select({}, "top5") == range(396, 400));
    CHECK(select({}, "top.500") == range(1, 400));
    CHECK(select({}, "excltop5") == range(1, 395));
    CHECK(select({}, "all") == range(1, 400));
    CHECK(select({}, "excl,all").empty());

    const LayerRules off = field("1-5,off");
    CHECK(off.off);
    CHECK(apply_rules(off).empty());
    const LayerRules nowarn = parse_layer_rules("-^ps=nowarn=ps^-", "nowarn.");
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
