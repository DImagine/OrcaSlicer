#ifndef slic3r_PreciseSeamRulesApply_hpp_
#define slic3r_PreciseSeamRulesApply_hpp_

#include <string>
#include <vector>
#include "libslic3r/Model.hpp"
#include "PreciseSeam.hpp"
#include "PreciseSeamRules.hpp"

// Precise Seam layer rules applied to a print: the branch that SeamPlacer and the G-code export take only
// for objects whose Precise Seam volumes use the rules language. Objects without rules keep the code
// path of Precise Seam without rules. Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class ConstLayerPtrsAdaptor;
namespace PreciseSeam {

// Whether any Precise Seam volume of the object uses the rules language: a non-empty rules field, or a
// name block marker in either direction, in any letter case. Broken or invalid rules count too, so that
// they are reported. The check reads the raw text and parses nothing.
bool uses_precise_seam_rules(const ModelObject &object);

// A rule diagnostic tied to its volume, for the export warnings.
struct ModifierDiagnostic {
    const ModelVolume *volume;
    RuleDiagnostic     diagnostic;
};

// init_precise_seam_data() for an object that uses layer rules: strong volumes in priority order, weak ones
// in application order, with their effective type and layers in `effective_out` (every Precise Seam volume
// that takes part, also those without rules of their own). Volumes switched off or excluded from every
// layer are left out. Diagnostics are appended to `diagnostics_out` in volume, source and element order.
// `model_layers` are the layers of the print object that `numbering` numbers, in the same order.
void init_precise_seam_data_by_rules(
    std::vector<const ModelVolume*> &strong_volumes_out,
    std::vector<const ModelVolume*> &weak_volumes_out,
    bool &has_strong_out,
    const ModelObject &model_object,
    const LayerNumbering &numbering,
    const ConstLayerPtrsAdaptor &model_layers,
    EffectiveModifiers &effective_out,
    std::vector<ModifierDiagnostic> &diagnostics_out);

// Single-line export warnings for the diagnostics of all print objects, in their order. A diagnostic that
// print objects of one model object repeat is reported once; equal texts of different volumes are merged
// with their count.
std::vector<std::string> rule_warnings(const std::vector<ModifierDiagnostic> &diagnostics);

} // namespace PreciseSeam
} // namespace Slic3r

#endif // slic3r_PreciseSeamRulesApply_hpp_
