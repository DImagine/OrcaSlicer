#include "PreciseSeamRules.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include "libslic3r/Model.hpp"
#include "libslic3r/libslic3r.h"

namespace Slic3r {
namespace PreciseSeam {

namespace {

constexpr std::string_view open_marker  = "{ps=";
constexpr std::string_view close_marker = "=ps}";

// ASCII only and independent of the locale, as the language is defined over ASCII.
char ascii_lower(char c) { return c >= 'A' && c <= 'Z' ? char(c - 'A' + 'a') : c; }
bool is_ascii_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\v' || c == '\f'; }

// Markers are matched without case so that names typed in CAD keep working.
bool marker_at(std::string_view text, size_t pos, std::string_view marker)
{
    if (text.size() - pos < marker.size())
        return false;
    for (size_t i = 0; i < marker.size(); ++i)
        if (ascii_lower(text[pos + i]) != marker[i])
            return false;
    return true;
}

struct NameBlocks {
    std::vector<std::string_view> blocks;
    std::vector<RuleDiagnostic>   diagnostics;
};

// Marker diagnostics refer to the name as a whole (name_block 0): a discarded block gets no number.
NameBlocks extract_name_blocks(std::string_view name)
{
    NameBlocks            out;
    std::optional<size_t> content_begin; // Set while inside a block.
    for (size_t pos = 0; pos < name.size();) {
        if (marker_at(name, pos, open_marker)) {
            // A new opening inside a block discards the started one.
            if (content_begin)
                out.diagnostics.push_back({RuleDiagnosticCode::BrokenMarker, {}, 0,
                                           std::string(name.substr(*content_begin, pos - *content_begin))});
            pos += open_marker.size();
            content_begin = pos;
        } else if (marker_at(name, pos, close_marker)) {
            if (content_begin) {
                out.blocks.push_back(name.substr(*content_begin, pos - *content_begin));
                content_begin.reset();
            } else
                out.diagnostics.push_back({RuleDiagnosticCode::BrokenMarker, {}, 0, {}}); // Stray closing marker.
            pos += close_marker.size();
        } else
            ++pos;
    }
    if (content_begin) // Unclosed at the end of the name.
        out.diagnostics.push_back({RuleDiagnosticCode::BrokenMarker, {}, 0, std::string(name.substr(*content_begin))});
    if (out.blocks.size() > 1)
        out.diagnostics.push_back({RuleDiagnosticCode::MultipleBlocks, {}, 0, {}});
    return out;
}

// Whitespace never separates elements: "1 2" is "12".
std::string normalize(std::string_view element)
{
    std::string out;
    out.reserve(element.size());
    for (char c : element)
        if (!is_ascii_space(c))
            out.push_back(ascii_lower(c));
    return out;
}

std::string_view trim(std::string_view text)
{
    while (!text.empty() && is_ascii_space(text.front()))
        text.remove_prefix(1);
    while (!text.empty() && is_ascii_space(text.back()))
        text.remove_suffix(1);
    return text;
}

// Reads one normalized element; a rule matches only if it consumes the element completely.
class Cursor
{
public:
    explicit Cursor(std::string_view text) : m_text(text) {}

    bool done() const { return m_pos == m_text.size(); }
    bool digit() const { return !done() && m_text[m_pos] >= '0' && m_text[m_pos] <= '9'; }
    bool out_of_range() const { return m_out_of_range; }

    bool eat(std::string_view word)
    {
        if (m_text.compare(m_pos, word.size(), word) != 0)
            return false;
        m_pos += word.size();
        return true;
    }
    // Takes the rest of the element as read.
    void skip_rest() { m_pos = m_text.size(); }
    bool eat(char c)
    {
        if (done() || m_text[m_pos] != c)
            return false;
        ++m_pos;
        return true;
    }
    // Leading zeros are allowed; the value saturates above the limit, which is remembered.
    std::optional<int64_t> uint()
    {
        const size_t begin = m_pos;
        int64_t      value = 0;
        while (!done() && m_text[m_pos] >= '0' && m_text[m_pos] <= '9') {
            value = std::min<int64_t>(value * 10 + (m_text[m_pos] - '0'), rule_number_limit + 1);
            ++m_pos;
        }
        if (m_pos == begin)
            return std::nullopt;
        if (value > rule_number_limit)
            m_out_of_range = true;
        return value;
    }
    // A relay-random seed, 0 .. 2^32 - 1, with leading zeros. Long input saturates; a larger value is no seed.
    std::optional<uint32_t> seed()
    {
        const size_t begin = m_pos;
        uint64_t     value = 0;
        for (; !done() && m_text[m_pos] >= '0' && m_text[m_pos] <= '9'; ++m_pos)
            value = std::min<uint64_t>(value * 10 + uint64_t(m_text[m_pos] - '0'), uint64_t(std::numeric_limits<uint32_t>::max()) + 1);
        if (m_pos == begin || value > std::numeric_limits<uint32_t>::max())
            return std::nullopt;
        return uint32_t(value);
    }
    // A decimal "uint [. uint]" in millionths: nanometres for a length in mm, millionths for a priority weight.
    // More than six decimals or a value above the limit is out of range, which is remembered; long input
    // saturates before that check.
    std::optional<int64_t> millionths()
    {
        const std::optional<int64_t> whole = uint();
        if (!whole)
            return std::nullopt;
        int64_t fraction = 0;
        if (eat('.')) {
            size_t digits = 0;
            for (; !done() && m_text[m_pos] >= '0' && m_text[m_pos] <= '9'; ++m_pos, ++digits)
                if (digits < 6)
                    fraction = fraction * 10 + (m_text[m_pos] - '0');
            if (digits == 0) // "12." is not a number.
                return std::nullopt;
            if (digits > 6)
                m_out_of_range = true;
            for (; digits < 6; ++digits)
                fraction *= 10;
        }
        // Six decimals: a millionth of a mm is a nanometre.
        static_assert(nanometres_per_mm == 1000000);
        const int64_t value = *whole * 1000000 + fraction;
        if (value > rule_number_limit * 1000000)
            m_out_of_range = true;
        return value;
    }

private:
    std::string_view m_text;
    size_t           m_pos          = 0;
    bool             m_out_of_range = false;
};

// Rounds towards positive infinity for a positive divisor.
int64_t ceil_div(int64_t numerator, int64_t divisor) { return numerator >= 0 ? (numerator + divisor - 1) / divisor : -(-numerator / divisor); }

// The rest of a height formula after "z" [sign] [k] "x": the levels k*x + b in mm for whole x >= 0, kept as the
// increasing levels above zero. Without such levels the set is valid and empty.
std::optional<RuleLayerSet> parse_height_formula(Cursor &c, int sign, std::optional<int64_t> lead)
{
    if (lead && *lead == 0) // "z0x+5" is a constant, not a formula.
        return std::nullopt;
    const int64_t k = (sign < 0 ? -1 : 1) * lead.value_or(nanometres_per_mm);
    int64_t       b = 0;
    const bool    plus = c.eat('+');
    if (plus || c.eat('-')) {
        const std::optional<int64_t> value = c.millionths();
        if (!value)
            return std::nullopt;
        b = plus ? *value : -*value;
    }
    // first > last: no level. With literals up to the limit, k*x0 stays within about 2e12.
    RuleLayerSet set{RuleLayerSet::Kind::Height, 1, 0, 1};
    if (k > 0) {
        // The smallest level above zero, without going back below b: x starts at 0.
        set.first = b + k * std::max<int64_t>(0, ceil_div(1 - b, k));
        set.last  = std::numeric_limits<int64_t>::max();
        set.step  = k;
    } else if (b > 0) {
        // b, b - |k|, ... down to the smallest level above zero.
        set.step                = -k;
        const int64_t remainder = b % set.step;
        set.first               = remainder != 0 ? remainder : set.step;
        set.last                = b;
    }
    return set;
}

// layer_set = uint | uint "-" uint [step] | formula | parity | top | all | height
std::optional<RuleLayerSet> parse_layer_set(Cursor &c)
{
    using Kind = RuleLayerSet::Kind;
    // height = "zin" decimal ["-" decimal] | "z" decimal ["-" decimal ["step" decimal]]
    //        | "z" ["+" | "-"] [decimal] "x" [("+" | "-") decimal]; zin is tried first, as z begins it.
    const bool inside = c.eat("zin");
    if (inside || c.eat('z')) {
        int sign = 0;
        if (!inside && c.eat('+'))
            sign = 1;
        else if (!inside && c.eat('-'))
            sign = -1;
        // A broken number is no omitted coefficient: "z1.x" is not "zx".
        std::optional<int64_t> first;
        if (c.digit() && !(first = c.millionths()))
            return std::nullopt;
        if (!inside && c.eat('x'))
            return parse_height_formula(c, sign, first);
        // Heights are unsigned: "z-1" and "z+1" are not valid.
        if (sign != 0 || !first)
            return std::nullopt;
        RuleLayerSet set{Kind::Height, *first, *first, 0};
        set.inside = inside;
        if (c.eat('-')) {
            const std::optional<int64_t> last = c.millionths();
            if (!last || *last < set.first)
                return std::nullopt;
            set.last = *last;
            if (!inside && c.eat("step")) {
                const std::optional<int64_t> step = c.millionths();
                if (!step || *step < 1)
                    return std::nullopt;
                set.step = *step;
            }
        }
        return set;
    }
    if (c.eat("odd"))
        return RuleLayerSet{Kind::Odd};
    if (c.eat("even"))
        return RuleLayerSet{Kind::Even};
    if (c.eat("all"))
        return RuleLayerSet{Kind::All};
    if (c.eat("top")) {
        const std::optional<int64_t> k = c.uint();
        if (!k || *k < 1)
            return std::nullopt;
        return RuleLayerSet{Kind::Top, *k};
    }

    // formula = [ "+" | "-" ], [ uint ], "x", [ ( "+" | "-" ), uint ]
    int sign = 0;
    if (c.eat('+'))
        sign = 1;
    else if (c.eat('-'))
        sign = -1;
    const std::optional<int64_t> lead = c.uint();
    if (c.eat('x')) {
        if (lead && *lead == 0) // "0x+5" is a constant, not a formula.
            return std::nullopt;
        RuleLayerSet set{Kind::Formula, (sign < 0 ? -1 : 1) * lead.value_or(1), 0};
        const bool plus = c.eat('+');
        if (plus || c.eat('-')) {
            const std::optional<int64_t> b = c.uint();
            if (!b)
                return std::nullopt;
            set.last = plus ? *b : -*b;
        }
        return set;
    }
    // Layer numbers and ranges are unsigned: "-5" and "+5" are not valid.
    if (sign != 0 || !lead || *lead < 1)
        return std::nullopt;
    RuleLayerSet set{Kind::Range, *lead, *lead};
    if (c.eat('-')) {
        const std::optional<int64_t> last = c.uint();
        if (!last || *last < set.first)
            return std::nullopt;
        set.last = *last;
        if (c.eat("step")) {
            const std::optional<int64_t> step = c.uint();
            if (!step || *step < 1)
                return std::nullopt;
            set.step = *step;
        }
    }
    return set;
}

// --- Word families -------------------------------------------------------------------------------------------
// Each family of words reads its own payload from the start of a normalized element and commits it to the rules.
// parse() returns NotMine when the element does not start with one of its keywords, and then reads nothing
// else; once its keyword matched, the element is the family's, valid or not. commit() runs only after the whole
// element was read without error, so a broken element changes nothing.

enum class ParseStatus { NotMine, Ok, Invalid };

// The state of one source while it is read: the mode and origin of the sets that follow. Each source starts
// in the include mode and counts from the bed.
struct SourceState {
    bool exclude = false;
    bool model   = false;
};

// One element as the families see it.
struct ElementRef {
    RuleSource       source;
    size_t           index = 0;     // Number in the source, from 1, empty elements counted.
    bool             first = false; // The first non-empty element of the source.
    std::string_view raw;           // Trimmed original text.
    std::string_view normalized;
};

void add_diagnostic(LayerRules &rules, RuleDiagnosticCode code, const ElementRef &element)
{
    rules.diagnostics.push_back({code, element.source, element.index, std::string(element.raw)});
}

void add_operation(LayerRules &rules, const SourceState &state, const RuleLayerSet &set, const ElementRef &element)
{
    rules.operations.push_back({state.exclude, state.model, set, element.source, element.index, std::string(element.raw)});
}

// The first value of a singleton word stays in effect; repeating an equal one is no conflict.
template<typename T>
void commit_first(std::optional<T> &slot, const T &value, RuleDiagnosticCode conflict, LayerRules &rules, const ElementRef &element)
{
    if (!slot)
        slot = value;
    else if (*slot != value)
        add_diagnostic(rules, conflict, element);
}

// id: a name for external scripts. Any text may follow, so it takes the element as a whole; no other keyword
// starts with "id". Repeats compare normalized, the trimmed original is kept.
struct IdWord {
    struct Payload {};
    static ParseStatus parse(Cursor &c, Payload &)
    {
        if (!c.eat("id"))
            return ParseStatus::NotMine;
        c.skip_rest();
        return ParseStatus::Ok;
    }
    static void commit(const Payload &, LayerRules &rules, SourceState &, const ElementRef &element)
    {
        if (!rules.id)
            rules.id = std::string(element.raw);
        else if (normalize(*rules.id) != element.normalized)
            add_diagnostic(rules, RuleDiagnosticCode::IdConflict, element);
    }
};

// comment = ";" {any character except ","}: a note for the reader, as in G-code; it changes nothing, may repeat
// and is never reported. No other keyword starts with ";".
struct CommentWord {
    struct Payload {};
    static ParseStatus parse(Cursor &c, Payload &)
    {
        if (!c.eat(';'))
            return ParseStatus::NotMine;
        c.skip_rest();
        return ParseStatus::Ok;
    }
    static void commit(const Payload &, LayerRules &, SourceState &, const ElementRef &) {}
};

// mode = ("incl" | "excl") [layer_set]: switches the mode of the source, and may add a set in it.
struct ModeWord {
    struct Payload {
        bool                        exclude = false;
        std::optional<RuleLayerSet> set;
    };
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        const bool incl = c.eat("incl");
        if (!incl && !c.eat("excl"))
            return ParseStatus::NotMine;
        out.exclude = !incl;
        if (!c.done() && !(out.set = parse_layer_set(c)))
            return ParseStatus::Invalid;
        return ParseStatus::Ok;
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &state, const ElementRef &element)
    {
        state.exclude = in.exclude;
        if (in.set)
            add_operation(rules, state, *in.set, element);
    }
};

// The rest of an apex after its keyword: [decimal], a threshold in mm.
ParseStatus parse_apex_rest(Cursor &c, ApexRule &out)
{
    if (!c.done() && !(out.threshold = c.millionths()))
        return ParseStatus::Invalid;
    return ParseStatus::Ok;
}

// type = label ["apex" [decimal]]: the first label stays in effect; gray and grey are one type. "center apex" is
// "center, apex" in one element; after another type the apex is reported and ignored once the type is known.
struct TypeWord {
    struct Payload {
        ModelVolumeType         type = ModelVolumeType::INVALID;
        std::optional<ApexRule> apex;
    };
    static constexpr std::pair<std::string_view, ModelVolumeType> labels[] = {
        {"center", ModelVolumeType::PRECISE_SEAM_CENTER},  {"left", ModelVolumeType::PRECISE_SEAM_LEFT},
        {"right", ModelVolumeType::PRECISE_SEAM_RIGHT},    {"green", ModelVolumeType::PRECISE_SEAM_ENFORCED},
        {"red", ModelVolumeType::PRECISE_SEAM_BLOCKED},    {"gray", ModelVolumeType::PRECISE_SEAM_NEUTRAL},
        {"grey", ModelVolumeType::PRECISE_SEAM_NEUTRAL},
    };
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        for (const auto &[label, type] : labels)
            if (c.eat(label)) {
                out.type = type;
                if (!c.eat("apex"))
                    return ParseStatus::Ok;
                return parse_apex_rest(c, out.apex.emplace());
            }
        return ParseStatus::NotMine;
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &, const ElementRef &element)
    {
        commit_first(rules.type, in.type, RuleDiagnosticCode::TypeConflict, rules, element);
        if (in.apex)
            commit_first(rules.apex, *in.apex, RuleDiagnosticCode::ApexConflict, rules, element);
    }
};

// priority = "pr" ["-neg"] decimal: a weight in millionths, the higher first; pr-negN weighs -N, and pr0, pr-neg0
// and no priority weigh the same. pr-N is no priority: a minus that only follows pr is easy to miss.
struct PriorityWord {
    struct Payload {
        int64_t priority = 0;
    };
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        if (!c.eat("pr"))
            return ParseStatus::NotMine;
        const bool                   negative = c.eat("-neg");
        const std::optional<int64_t> value    = c.millionths();
        if (!value)
            return ParseStatus::Invalid;
        out.priority = negative ? -*value : *value;
        return ParseStatus::Ok;
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &, const ElementRef &element)
    {
        commit_first(rules.priority, in.priority, RuleDiagnosticCode::PriorityConflict, rules, element);
    }
};

// Flags of the volume, which may repeat: off, blank, nowarn.
struct FlagWord {
    struct Payload {
        bool LayerRules::*flag = nullptr;
    };
    static constexpr std::pair<std::string_view, bool LayerRules::*> flags[] = {
        {"off", &LayerRules::off}, {"blank", &LayerRules::blank}, {"nowarn", &LayerRules::nowarn}};
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        for (const auto &[word, flag] : flags)
            if (c.eat(word)) {
                out.flag = flag;
                return ParseStatus::Ok;
            }
        return ParseStatus::NotMine;
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &, const ElementRef &) { rules.*in.flag = true; }
};

// override: only as the first element of the field, where parse_layer_rules() already skipped the name.
struct OverrideWord {
    struct Payload {};
    static ParseStatus parse(Cursor &c, Payload &) { return c.eat("override") ? ParseStatus::Ok : ParseStatus::NotMine; }
    // The whole normalized element is override, as read_word() would accept it; for the look ahead at the name.
    static bool whole(std::string_view normalized)
    {
        Cursor  c(normalized);
        Payload payload;
        return parse(c, payload) == ParseStatus::Ok && c.done() && !c.out_of_range();
    }
    static void commit(const Payload &, LayerRules &rules, SourceState &, const ElementRef &element)
    {
        if (element.source.field && element.first)
            rules.name_overridden = true;
        else
            add_diagnostic(rules, RuleDiagnosticCode::InvalidElement, element);
    }
};

// relay = "relay" [uint] | "relay-random" [uint] ["seed" seed]; without a number it is a group of its own,
// apart from relay1, ... relay-random is tried first, as relay begins it.
struct RelayWord {
    struct Payload {
        RelayLabel              label;
        std::optional<uint32_t> seed;
    };
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        const bool random = c.eat("relay-random");
        if (!random && !c.eat("relay"))
            return ParseStatus::NotMine;
        out.label.random = random;
        bool seed        = random && c.eat("seed");
        if (!seed && !c.done()) {
            const std::optional<int64_t> group = c.uint();
            if (!group || *group < 1)
                return ParseStatus::Invalid;
            out.label.number = *group;
            seed             = random && c.eat("seed");
        }
        if (seed && !(out.seed = c.seed()))
            return ParseStatus::Invalid;
        return ParseStatus::Ok;
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &, const ElementRef &element)
    {
        // The seed is a part of the element: a later relay with another seed, or without one, conflicts.
        if (!rules.relay) {
            rules.relay      = in.label;
            rules.relay_seed = in.seed;
        } else if (*rules.relay != in.label || rules.relay_seed != in.seed)
            add_diagnostic(rules, RuleDiagnosticCode::RelayConflict, element);
    }
};

// origin = "model" | "bed": the sets after it in the source count model layers and heights, or those of the bed.
struct OriginWord {
    struct Payload {
        bool model = false;
    };
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        out.model = c.eat("model");
        return out.model || c.eat("bed") ? ParseStatus::Ok : ParseStatus::NotMine;
    }
    static void commit(const Payload &in, LayerRules &, SourceState &state, const ElementRef &) { state.model = in.model; }
};

// magnet = "magnet" [decimal] ["apex" [decimal]]; distances in mm are kept in nanometres, 0 means a touch.
struct MagnetWord {
    using Payload = MagnetRule;
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        if (!c.eat("magnet"))
            return ParseStatus::NotMine;
        out.apex = c.eat("apex");
        if (!out.apex && !c.done()) {
            if (!(out.limit = c.millionths()))
                return ParseStatus::Invalid;
            out.apex = c.eat("apex");
        }
        if (out.apex && !c.done()) {
            const std::optional<int64_t> length = c.millionths();
            if (!length)
                return ParseStatus::Invalid;
            out.apex_length = *length;
        }
        return ParseStatus::Ok;
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &, const ElementRef &element)
    {
        // Distances compare as numbers: magnet2.5 repeats magnet02.500000; magnet2apex differs from magnet2.
        commit_first(rules.magnet, in, RuleDiagnosticCode::MagnetConflict, rules, element);
    }
};

// apex = "apex" [decimal]: the same as after a type; apart from "magnet ... apex", which belongs to the magnet.
struct ApexWord {
    using Payload = ApexRule;
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        if (!c.eat("apex"))
            return ParseStatus::NotMine;
        return parse_apex_rest(c, out);
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &, const ElementRef &element)
    {
        // Thresholds compare as numbers: apex1 repeats apex1.0; apex and apex0.012 differ.
        commit_first(rules.apex, in, RuleDiagnosticCode::ApexConflict, rules, element);
    }
};

// A layer set on its own, added in the current mode. Last: whatever no keyword claimed is read as a set.
struct SetWord {
    using Payload = RuleLayerSet;
    static ParseStatus parse(Cursor &c, Payload &out)
    {
        std::optional<RuleLayerSet> set = parse_layer_set(c);
        if (!set)
            return ParseStatus::Invalid;
        out = std::move(*set);
        return ParseStatus::Ok;
    }
    static void commit(const Payload &in, LayerRules &rules, SourceState &state, const ElementRef &element)
    {
        add_operation(rules, state, in, element);
    }
};

// Reads the element with one family on a fresh cursor and commits it when it is the family's and valid.
template<typename Word>
ParseStatus read_word(LayerRules &rules, SourceState &state, const ElementRef &element)
{
    Cursor                 c(element.normalized);
    typename Word::Payload payload{};
    const ParseStatus      status = Word::parse(c, payload);
    if (status == ParseStatus::NotMine)
        return status;
    // A malformed element changes nothing, not even the mode.
    if (status == ParseStatus::Invalid || !c.done() || c.out_of_range()) {
        add_diagnostic(rules, c.out_of_range() ? RuleDiagnosticCode::NumberOutOfRange : RuleDiagnosticCode::InvalidElement, element);
        return ParseStatus::Invalid;
    }
    Word::commit(payload, rules, state, element);
    return ParseStatus::Ok;
}

// The families in the order they are tried. No keyword of the language is a prefix of another but those noted
// at their families; a comment and id come first, as any text may follow them, and a bare layer set last.
using WordReader = ParseStatus (*)(LayerRules &, SourceState &, const ElementRef &);
constexpr WordReader words[] = {read_word<CommentWord>,  read_word<IdWord>,
                                read_word<ModeWord>,     read_word<TypeWord>,
                                read_word<PriorityWord>, read_word<FlagWord>,  read_word<OverrideWord>,
                                read_word<RelayWord>,    read_word<OriginWord>, read_word<MagnetWord>,
                                read_word<ApexWord>,     read_word<SetWord>};

// Reads one source: sets, type and priority carry over between sources, the mode and origin do not.
void parse_source(LayerRules &rules, std::string_view source_text, const RuleSource &source)
{
    SourceState state;
    bool        first = true; // No non-empty element seen yet.
    size_t      index = 0;
    for (size_t begin = 0; begin <= source_text.size(); ++index) {
        const size_t           comma = std::min(source_text.find(',', begin), source_text.size());
        const std::string_view raw   = trim(source_text.substr(begin, comma - begin));
        begin                        = comma + 1;

        const std::string normalized = normalize(raw);
        if (normalized.empty())
            continue;
        const ElementRef element{source, index + 1, std::exchange(first, false), raw, normalized};
        // The set family comes last and claims every element, so one family always decides.
        for (const WordReader read : words)
            if (read(rules, state, element) != ParseStatus::NotMine)
                break;
    }
}

// The first element that is not empty after normalization, as parse_source() reads it.
std::string leading_element(std::string_view source_text)
{
    for (size_t begin = 0; begin <= source_text.size();) {
        const size_t comma      = std::min(source_text.find(',', begin), source_text.size());
        std::string  normalized = normalize(source_text.substr(begin, comma - begin));
        if (!normalized.empty())
            return normalized;
        begin = comma + 1;
    }
    return {};
}

// Calls `visit` for every layer of the set within [1, count]. With literals up to rule_number_limit, every
// layer value of the numbered sets stays within +-(count + 3 * rule_number_limit): a formula starts at b or
// within one step k of 1 or count, and every loop stops at most one step past its bound. Height levels have
// their own bounds, below.
template<class Visit>
void for_each_layer(const RuleLayerSet &set, int64_t count, const std::vector<size_t> &object_layers, const std::vector<LayerBand> &bands,
                    Visit &&visit)
{
    using Kind = RuleLayerSet::Kind;
    switch (set.kind) {
    case Kind::Range:
        for (int64_t layer = set.first; layer <= std::min(set.last, count); layer += set.step)
            visit(layer);
        break;
    case Kind::Odd:
    case Kind::Even:
        for (int64_t layer = set.kind == Kind::Odd ? 1 : 2; layer <= count; layer += 2)
            visit(layer);
        break;
    case Kind::All:
        for (int64_t layer = 1; layer <= count; ++layer)
            visit(layer);
        break;
    case Kind::Top: {
        // The highest numbered layers of this object, regardless of taller objects on the plate.
        auto end = std::upper_bound(object_layers.begin(), object_layers.end(), size_t(count));
        for (auto it = end - std::min<int64_t>(set.first, end - object_layers.begin()); it != end; ++it)
            if (*it >= 1)
                visit(int64_t(*it));
        break;
    }
    case Kind::Formula: {
        // Only x >= 0 whose value k*x+b falls into [1, count] is visited, without scanning from x = 0.
        const int64_t k = set.first, b = set.last;
        if (k > 0) {
            for (int64_t layer = k * std::max<int64_t>(0, ceil_div(1 - b, k)) + b; layer <= count; layer += k)
                visit(layer);
        } else {
            for (int64_t layer = b - (-k) * std::max<int64_t>(0, ceil_div(b - count, -k)); layer >= 1; layer += k)
                visit(layer);
        }
        break;
    }
    case Kind::Height: {
        if (set.step > 0) {
            if (set.first > set.last)
                break;
            // Levels first + i*step, i in [0, n], computed exactly; first >= 0, so n cannot overflow. Each
            // level selects as zA does. The first level above a band's bottom is found by a binary search, as
            // its predicate is monotonic in i also after rounding; when it does not reach into the band, no
            // later level does. This takes no time per level, whatever their number.
            const int64_t n     = (set.last - set.first) / set.step;
            const auto    level = [&set](int64_t i) { return double(set.first + i * set.step) / nanometres_per_mm; };
            for (const LayerBand &band : bands) {
                if (band.number < 1 || band.number > size_t(count) || !(band.bottom_z < level(n) - EPSILON))
                    continue;
                int64_t lo = 0, hi = n;
                while (lo < hi) {
                    const int64_t mid = lo + (hi - lo) / 2;
                    if (band.bottom_z < level(mid) - EPSILON)
                        hi = mid;
                    else
                        lo = mid + 1;
                }
                if (band.print_z >= level(lo) - EPSILON)
                    visit(int64_t(band.number));
            }
            break;
        }
        // A whole layer never fits into a single height, whatever the tolerance.
        if (set.inside && set.first == set.last)
            break;
        // Bands are (bottom_z, print_z]. Both bounds of an intersection lean towards the lower layer by EPSILON.
        const double a = double(set.first) / nanometres_per_mm, b = double(set.last) / nanometres_per_mm;
        for (const LayerBand &band : bands)
            if (band.number >= 1 && band.number <= size_t(count) &&
                (set.inside ? band.bottom_z >= a - EPSILON && band.print_z <= b + EPSILON :
                              band.print_z >= a - EPSILON && band.bottom_z < b - EPSILON))
                visit(int64_t(band.number));
        break;
    }
    }
}

} // namespace

std::string relay_label_text(const RelayLabel &label)
{
    std::string text = label.random ? "relay-random" : "relay";
    if (label.number != 0)
        text += std::to_string(label.number);
    return text;
}

bool has_rule_marker(std::string_view volume_name)
{
    for (size_t pos = 0; pos < volume_name.size(); ++pos)
        if (marker_at(volume_name, pos, open_marker) || marker_at(volume_name, pos, close_marker))
            return true;
    return false;
}

bool uses_precise_seam_rules(const ModelObject &object)
{
    return std::any_of(object.volumes.begin(), object.volumes.end(), [](const ModelVolume *volume) {
        return volume->is_precise_seam() && (!volume->precise_seam_rules.empty() || has_rule_marker(volume->name));
    });
}

LayerRules parse_layer_rules(std::string_view volume_name, std::string_view field)
{
    LayerRules rules;
    // override first in the field replaces the name: neither its rules nor its diagnostics take part.
    if (!OverrideWord::whole(leading_element(field))) {
        NameBlocks name   = extract_name_blocks(volume_name);
        rules.diagnostics = std::move(name.diagnostics);
        for (size_t i = 0; i < name.blocks.size(); ++i)
            parse_source(rules, name.blocks[i], {false, i + 1});
    }
    parse_source(rules, field, {true, 0});
    return rules;
}

LayerSelection select_layers(const LayerRules &rules, size_t layer_count, const std::vector<size_t> &object_layers,
                             const std::vector<LayerBand> &bands, const std::vector<size_t> &model_to_bed, double model_bottom)
{
    LayerSelection selection;
    selection.layers.assign(layer_count + 1, false);
    if (rules.inert())
        return selection;
    // No run-time overflow check, by design: a check is kept only for a case that can happen. layer_count is
    // the size of the vector just allocated, so the arithmetic of for_each_layer() stays far inside int64_t;
    // the assert states this invariant in debug builds.
    assert(layer_count <= size_t(std::numeric_limits<int64_t>::max() - 3 * rule_number_limit));
    const int64_t count = int64_t(layer_count);
    // Without any valid set the modifier keeps its behaviour on every layer.
    // Otherwise the mode of the first valid set picks the initial set: all layers for exclusion.
    if (!rules.restricts_layers() || rules.operations.front().exclude)
        std::fill(selection.layers.begin() + 1, selection.layers.begin() + 1 + count, true);
    // Model heights: the same bands, measured from the bottom of the model.
    std::vector<LayerBand> model_bands;
    for (const RuleLayerOperation &operation : rules.operations)
        if (operation.model && operation.set.kind == RuleLayerSet::Kind::Height) {
            for (const LayerBand &band : bands)
                model_bands.push_back({band.number, band.bottom_z - model_bottom, band.print_z - model_bottom});
            break;
        }
    for (const RuleLayerOperation &operation : rules.operations) {
        bool       any   = false;
        const auto visit = [&](int64_t layer) {
            selection.layers[size_t(layer)] = !operation.exclude;
            any                             = true;
        };
        using Kind = RuleLayerSet::Kind;
        if (!operation.model || operation.set.kind == Kind::Top)
            // topK takes the highest numbered model layers in either origin.
            for_each_layer(operation.set, count, object_layers, bands, visit);
        else if (operation.set.kind == Kind::Height)
            for_each_layer(operation.set, count, object_layers, model_bands, visit);
        else
            // Model layer i is the number of the i-th model layer; one without a number selects nothing, and
            // the numbers of the model layers above it do not move down.
            for_each_layer(operation.set, int64_t(model_to_bed.size()), {}, {}, [&](int64_t layer) {
                const size_t number = model_to_bed[size_t(layer) - 1];
                assert(number <= layer_count);
                if (number != 0)
                    visit(int64_t(number));
            });
        if (!any)
            selection.diagnostics.push_back({RuleDiagnosticCode::EmptySet, operation.source, operation.element, operation.text});
    }
    return selection;
}

std::vector<size_t> LayerNumbering::numbered_layers() const
{
    std::vector<size_t> out;
    out.reserve(object_layers.size());
    for (size_t number : object_layers)
        if (number != 0)
            out.push_back(number);
    // Numbers grow with the object layers; sorting only guards the contract of select_layers().
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

EffectiveOrder effective_order(const std::vector<RuledModifier> &modifiers)
{
    EffectiveOrder order;
    for (size_t i = 0; i < modifiers.size(); ++i) {
        const RuledModifier &modifier = modifiers[i];
        if (modifier.rules->inert())
            continue;
        const ModelVolumeType type = modifier.effective_type();
        (is_precise_seam_strong(type) ? order.strong : order.weak).push_back(i);
        if (type != modifier.actual)
            (modifier.rules->nowarn ? order.suppressed : order.type_changed).push_back(i);
    }
    // The higher weight first; a volume without priority weighs 0, as pr0.
    const auto weight = [&modifiers](size_t i) { return modifiers[i].rules->priority.value_or(0); };
    std::vector<bool> reordered(modifiers.size(), false);
    for (std::vector<size_t> *group : {&order.strong, &order.weak}) {
        // Without priorities the tree order stands and nothing is reordered.
        if (std::none_of(group->begin(), group->end(), [&modifiers](size_t i) { return modifiers[i].rules->priority.has_value(); }))
            continue;
        // The group is collected in input order, so a stable sort breaks ties by the input index.
        const std::vector<size_t> tree = *group;
        std::stable_sort(group->begin(), group->end(), [&weight](size_t a, size_t b) { return weight(a) > weight(b); });
        // A pair whose order differs from the tree reports its members that carry a priority.
        std::vector<size_t> position(modifiers.size());
        for (size_t p = 0; p < group->size(); ++p)
            position[(*group)[p]] = p;
        for (size_t a = 0; a < tree.size(); ++a)
            for (size_t b = a + 1; b < tree.size(); ++b)
                if (position[tree[b]] < position[tree[a]])
                    for (size_t i : {tree[a], tree[b]})
                        if (modifiers[i].rules->priority)
                            reordered[i] = true;
    }
    for (size_t i = 0; i < modifiers.size(); ++i)
        if (reordered[i]) {
            if (!modifiers[i].rules->nowarn)
                order.reordered.push_back(i);
            else if (std::find(order.suppressed.begin(), order.suppressed.end(), i) == order.suppressed.end())
                order.suppressed.push_back(i);
        }
    std::sort(order.suppressed.begin(), order.suppressed.end());
    return order;
}

} // namespace PreciseSeam
} // namespace Slic3r
