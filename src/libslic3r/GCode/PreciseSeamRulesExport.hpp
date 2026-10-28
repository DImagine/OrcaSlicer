#ifndef slic3r_PreciseSeamRulesExport_hpp_
#define slic3r_PreciseSeamRulesExport_hpp_

#include <functional>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>
#include "libslic3r/GCode.hpp"
#include "libslic3r/libslic3r.h"
#include "PreciseSeamRules.hpp"
#include "SeamPlacer.hpp"

// Precise Seam layer rules in the G-code export: the layer plan collected early, once, and the rule numbers of
// the layers taken from it. Only GCode.cpp includes this header, so that the rules headers that SeamPlacer.hpp
// reaches never include GCode.hpp. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class Layer;
class Print;
class PrintObject;
class ToolOrdering;
namespace PreciseSeam {

class ExportPlan
{
public:
    using LayersByZ    = std::vector<std::pair<coordf_t, std::vector<GCode::LayerToPrint>>>;
    using ObjectLayers = std::vector<GCode::LayerToPrint>;

    // What the export lends for the preparation, which only uses them before it returns.
    struct Access {
        std::function<LayersByZ()>                              collect_print;        // GCode::collect_layers_to_print(print).
        std::function<ObjectLayers(const PrintObject &)>        collect_object;       // Of one object, with the belt flag.
        std::function<bool(const PrintObject &, const Layer &)> object_layer_prints;  // A belt object layer writes something.
    };

    // Collects the plan of the export, only when some object uses rules, and numbers its layers. `by_object`
    // is the export's own condition for printing object by object: the by-object sequence without a wipe tower.
    void prepare(const Print &print, const ToolOrdering &tool_ordering, bool by_object, const Access &access);

    // Rule numbers of the layers of each object that uses rules, for SeamPlacer::init().
    const SeamPlacer::PreciseSeamNumbering &numbering() const { return m_numbering; }
    // The plan by layer collected here, moved out, or nullopt when there is none: collect_layers_to_print() raises
    // warnings and errors, so the export collects it only then.
    std::optional<LayersByZ> take_layers_by_z() { return std::exchange(m_by_z, std::nullopt); }
    // The plan of one object collected here, or nullptr when there is none.
    const ObjectLayers *object_layers(const PrintObject &object) const;

private:
    std::optional<LayersByZ>                                m_by_z;
    std::unordered_map<const PrintObject *, ObjectLayers>   m_by_object;
    SeamPlacer::PreciseSeamNumbering                        m_numbering;
};

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_PreciseSeamRulesExport_hpp_
