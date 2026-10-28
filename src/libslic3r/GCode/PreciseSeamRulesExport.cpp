// Precise Seam layer rules in the G-code export. Design: docs/HLSD/precise-seam.md
#include "PreciseSeamRulesExport.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include "libslic3r/BeltBrim.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"
#include "PreciseSeamRules.hpp"
#include "PreciseSeamRulesApply.hpp"
#include "ToolOrdering.hpp"

namespace Slic3r {
namespace PreciseSeam {

namespace {

// Precise Seam layer rules count the layers as the preview slider does: one number per event of the
// export plan that prints something. An event prints when its tools have extruders, as process_layer()
// requires, or when it is a belt brim apron layer, which is handled before that check.
SeamPlacer::PreciseSeamNumbering print_numbering(const Print &print, const ToolOrdering &tool_ordering,
                                                 const ExportPlan::LayersByZ &layers_to_print)
{
    // Rows only for the objects with rules; the events of all objects take numbers.
    SeamPlacer::PreciseSeamNumbering numbering;
    for (const PrintObject *object : print.objects())
        if (uses_precise_seam_rules(*object->model_object()))
            numbering[object].object_layers.assign(object->layers().size(), 0);
    size_t count = 0;
    for (const auto &[print_z, layers] : layers_to_print) {
        const bool belt_brim_only =
            std::none_of(layers.begin(), layers.end(), [](const GCode::LayerToPrint &l) { return l.object_layer || l.support_layer; }) &&
            std::any_of(layers.begin(), layers.end(), [](const GCode::LayerToPrint &l) { return l.belt_brim_band != nullptr; });
        if (!belt_brim_only && tool_ordering.tools_for_layer(print_z).extruders.empty())
            continue;
        ++count;
        // Entries follow print.objects(); an event can hold the support of one object and a model layer of another.
        for (size_t i = 0; i < layers.size(); ++i)
            if (const Layer *layer = layers[i].object_layer)
                if (const auto row = numbering.find(print.objects()[i]); row != numbering.end()) {
                    const size_t index = layer->id() - row->first->slicing_parameters().raft_layers();
                    if (index < row->second.object_layers.size())
                        row->second.object_layers[index] = count;
                }
    }
    for (auto &entry : numbering)
        entry.second.count = count;
    return numbering;
}

// Printing by object: the object's own events from 1, the same for all its copies, so that they share
// one seam cache. Only events that write something are numbered, as when printing by layer.
LayerNumbering own_numbering(const PrintObject &object, const ExportPlan::ObjectLayers &layers_to_print,
                                const std::function<bool(const PrintObject &, const Layer &)> &object_layer_prints)
{
    LayerNumbering numbering;
    numbering.object_layers.assign(object.layers().size(), 0);
    for (const GCode::LayerToPrint &ltp : layers_to_print) {
        const bool prints = (ltp.object_layer && object_layer_prints(object, *ltp.object_layer)) ||
                            (ltp.support_layer && ltp.support_layer->has_extrusions()) ||
                            (ltp.belt_brim_band && !ltp.belt_brim_band->fills.empty());
        if (!prints)
            continue;
        ++numbering.count;
        if (ltp.object_layer) {
            const size_t index = ltp.object_layer->id() - object.slicing_parameters().raft_layers();
            if (index < numbering.object_layers.size())
                numbering.object_layers[index] = numbering.count;
        }
    }
    return numbering;
}

} // namespace

void ExportPlan::prepare(const Print &print, const ToolOrdering &tool_ordering, bool by_object, const Access &access)
{
    const auto uses_rules = [](const PrintObject *object) { return uses_precise_seam_rules(*object->model_object()); };
    if (by_object) {
        for (const PrintObject *object : print.objects())
            if (uses_rules(object)) {
                const ObjectLayers &layers = m_by_object[object] = access.collect_object(*object);
                m_numbering[object]        = own_numbering(*object, layers, access.object_layer_prints);
            }
    } else if (std::any_of(print.objects().begin(), print.objects().end(), uses_rules)) {
        // Printing by layer, the numbers count the events of the whole plate, objects without rules included.
        m_by_z      = access.collect_print();
        m_numbering = print_numbering(print, tool_ordering, *m_by_z);
    }
}

const ExportPlan::ObjectLayers *ExportPlan::object_layers(const PrintObject &object) const
{
    const auto it = m_by_object.find(&object);
    return it != m_by_object.end() ? &it->second : nullptr;
}

} // namespace PreciseSeam
} // namespace Slic3r
