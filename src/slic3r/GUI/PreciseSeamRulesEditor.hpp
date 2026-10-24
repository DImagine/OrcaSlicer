#ifndef slic3r_PreciseSeamRulesEditor_hpp_
#define slic3r_PreciseSeamRulesEditor_hpp_

#include <string>
#include <vector>

#include <wx/weakref.h>

#include "libslic3r/ObjectID.hpp"

class wxSizer;
class wxTextCtrl;
class wxWindow;

// Precise Seam rules field of the part settings: a copy on the Frequent and the Others page of the part
// tab, shown for a single selected Precise Seam volume. The volume holds the only accepted text.
// Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class ModelVolume;
namespace GUI {

// Writes raw rules to the volume with `volume_id` in the current model: when the text differs, one undo
// snapshot, then a background update. False when that volume is gone or no longer Precise Seam.
bool apply_precise_seam_rules(ObjectID volume_id, const std::string &rules);

class PreciseSeamRulesEditor
{
public:
    // The volume whose rules the fields show; nullptr hides them. A pending edit goes to its own volume first.
    void set_volume(const ModelVolume *volume);
    bool shown() const { return m_volume.valid(); }

    // Builds one copy of the field; pages build them lazily and destroy them when another page is shown.
    wxSizer *create(wxWindow *parent);

    // Applies the pending edit to the volume it was typed for, whatever is selected now.
    void commit_pending();
    // Drops the pending edit and shows the model, after undo, redo or a new project replaced it.
    void reset_from_model();
    // Shows the accepted text, or the pending edit of the shown volume.
    void refresh();

private:
    void show(const std::string &text);

    ObjectID                         m_volume;  // Volume shown; invalid when hidden.
    ObjectID                         m_edited;  // Volume of the pending edit; invalid without one.
    std::string                      m_draft;   // Text of the pending edit.
    std::vector<wxWeakRef<wxTextCtrl>> m_fields; // Live copies; a destroyed one turns null.
    bool                             m_showing = false; // Programmatic update in progress.
};

// The editor of the part settings tab, for callers outside the tab; no-ops before the tab exists.
void commit_precise_seam_rules_edit();
void reset_precise_seam_rules_edit();

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_PreciseSeamRulesEditor_hpp_
