#ifndef slic3r_PreciseSeamRulesEditor_hpp_
#define slic3r_PreciseSeamRulesEditor_hpp_

#include <string>
#include <vector>

#include <wx/string.h>
#include <wx/weakref.h>

#include "libslic3r/ObjectID.hpp"

class wxSizer;
class wxTextCtrl;
class wxWindow;

// Precise Seam rules field of the part settings: a copy on the Frequent and the Others page of the part
// tab, shown for a single selected Precise Seam volume. Besides the volume name, it is the only editor of
// the rules; the volume holds the only accepted text.
// Design: docs/HLSD/precise-seam.md

namespace Slic3r {
class ModelVolume;
namespace GUI {
class Page;

// Help of the rules language, shown on the title of the settings group.
wxString precise_seam_rules_help();

// Writes raw rules to the volume with `volume_id` in the current model: when the text differs, one undo
// snapshot, then a background update. False when that volume is gone or no longer Precise Seam.
bool apply_precise_seam_rules(ObjectID volume_id, const std::string &rules);

class PreciseSeamRulesEditor
{
public:
    // Adds a settings group with one copy of the field to `page`; the group shows while a volume is shown.
    void add_group(Page &page);
    // Shows the rules of the selected part when it is a single Precise Seam volume, otherwise hides them;
    // `single` is null for any other selection. A pending edit goes to its own volume first. True when the
    // fields appear or disappear, so the page must update its visibility.
    bool select(const ObjectBase *single);
    bool shown() const { return m_volume.valid(); }
    // Typed text not applied yet; saving and closing must see it as a change.
    bool pending() const { return m_edited.valid(); }

    // Builds one copy of the field; pages build them lazily and destroy them when another page is shown.
    wxSizer *create(wxWindow *parent);

    // Applies the pending edit to the volume it was typed for, whatever is selected now.
    void commit_pending();
    // Drops the pending edit and shows the model, after undo, redo or a new project replaced it.
    void reset_from_model();
    // Shows the accepted text, or the pending edit of the shown volume.
    void refresh();
    // Height and colours of the live copies after a DPI or a theme change.
    void rescale();

private:
    void show(const std::string &text);
    void fit(wxTextCtrl *text) const;

    ObjectID                         m_volume;  // Volume shown; invalid when hidden.
    ObjectID                         m_edited;  // Volume of the pending edit; invalid without one.
    std::string                      m_draft;   // Text of the pending edit.
    std::vector<wxWeakRef<wxTextCtrl>> m_fields; // Live copies; a destroyed one turns null.
    bool                             m_showing = false; // Programmatic update in progress.
};

// The editor of the part settings tab, for callers outside the tab; no-ops before the tab exists.
void commit_precise_seam_rules_edit();
void reset_precise_seam_rules_edit();
bool precise_seam_rules_edit_pending();

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_PreciseSeamRulesEditor_hpp_
