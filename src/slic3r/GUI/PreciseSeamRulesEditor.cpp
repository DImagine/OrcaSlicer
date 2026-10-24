#include "PreciseSeamRulesEditor.hpp"

#include <algorithm>

#include <wx/sizer.h>
#include <wx/textctrl.h>

#include "libslic3r/Model.hpp"
#include "GUI.hpp"
#include "GUI_App.hpp"
#include "I18N.hpp"
#include "Plater.hpp"
#include "Tab.hpp"

namespace Slic3r {
namespace GUI {

namespace {
ModelVolume *find_volume(ObjectID volume_id)
{
    for (ModelObject *object : wxGetApp().model().objects)
        for (ModelVolume *volume : object->volumes)
            if (volume->id() == volume_id)
                return volume;
    return nullptr;
}

PreciseSeamRulesEditor *part_editor()
{
    // The project is reset before the settings tabs exist and after they are gone.
    if (wxGetApp().model_tabs_list.size() < 2)
        return nullptr;
    auto *tab = dynamic_cast<TabPrintPart *>(wxGetApp().get_model_tab(true));
    return tab ? &tab->precise_seam_rules() : nullptr;
}
} // namespace

wxString precise_seam_rules_help()
{
    // TRN: Help for the Precise Seam rules field. Keep the rule keywords untranslated.
    return _L("Experimental: the rules syntax may change in future versions.\n"
              "Comma-separated rules, applied from left to right:\n"
              "layers: 5, 10-20, 10-20step3, 3x+1, odd, even, top5, all\n"
              "heights in mm: z12.5-30 (layers that reach into it), zin12.5-30 (layers wholly inside), z12.5; "
              "z0-30step2.5 and z2.5x+0.4 select the layer at each level, as z12.5 does\n"
              "model, bed: the layers and heights after it count from the model, above a raft, or from the bed (default)\n"
              "mode: incl (default), excl\n"
              "type: center, left, right, green, red, gray\n"
              "priority: pr1 (highest), pr-1 (lowest); off; nowarn\n"
              "relay, relay1: modifiers of one group act in turn, layer by layer: all center, left, right, or all green, red, gray\n"
              "relay-random, relay-random1: the same, but each layer goes to a random member, not the previous one when another can; "
              "relay-random1 seed 7 takes another draw, which may still give the same members\n"
              "magnet, magnet2.5: where the modifier does not cross a wall, it acts on the point nearest to it; within 2.5 mm if given, magnet0 on a touch\n"
              "idName: a name for scripts; the rules ignore it\n"
              "Layers are numbered as on the preview slider; when printing by object, from 1 for each object.\n"
              "Heights are measured from the bed by default, as on the slider. Decimals take a point: a comma separates rules.\n"
              "The name may hold more rules between -^ps= and =ps^-; override as the first rule here ignores them.");
}

bool apply_precise_seam_rules(ObjectID volume_id, const std::string &rules)
{
    ModelVolume *volume = find_volume(volume_id);
    if (volume == nullptr || !volume->is_precise_seam())
        return false;
    if (rules != volume->precise_seam_rules) {
        wxGetApp().plater()->take_snapshot(_u8L("Edit Precise Seam Rules"));
        volume->precise_seam_rules = rules;
        // Rules change only the seam placement; the export step is invalidated when the print is updated.
        wxGetApp().plater()->schedule_background_process();
    }
    if (PreciseSeamRulesEditor *editor = part_editor())
        editor->refresh();
    return true;
}

void PreciseSeamRulesEditor::set_volume(const ModelVolume *volume)
{
    const ObjectID id = volume ? volume->id() : ObjectID();
    // The same volume may hold other rules now, after undo or redo; refresh() keeps a pending edit.
    if (id == m_volume) {
        refresh();
        return;
    }
    commit_pending();
    m_volume = id;
    refresh();
}

wxSizer *PreciseSeamRulesEditor::create(wxWindow *parent)
{
    auto *text = new wxTextCtrl(parent, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE);
    fit(text);
    text->SetToolTip(_L("Precise Seam rules of this volume; the group title shows the help. "
                        "They apply when the field loses focus or on Ctrl+Enter."));

    // ChangeValue() sends no text event, so only typing reaches here.
    text->Bind(wxEVT_TEXT, [this, text](wxCommandEvent &) {
        if (m_showing || !m_volume.valid())
            return;
        m_edited = m_volume;
        m_draft  = into_u8(text->GetValue());
        // The other copy, if it is built, shows the same draft.
        for (wxTextCtrl *field : m_fields)
            if (field != nullptr && field != text && field->GetValue() != text->GetValue())
                field->ChangeValue(text->GetValue());
    });
    text->Bind(wxEVT_KILL_FOCUS, [this](wxFocusEvent &event) {
        commit_pending();
        event.Skip();
    });
    text->Bind(wxEVT_KEY_DOWN, [this](wxKeyEvent &event) {
        // Enter alone starts a new line; Ctrl+Enter (Cmd+Enter on macOS) applies.
        if ((event.GetKeyCode() == WXK_RETURN || event.GetKeyCode() == WXK_NUMPAD_ENTER) && event.GetModifiers() == wxMOD_CONTROL)
            commit_pending();
        else
            event.Skip();
    });

    m_fields.erase(std::remove_if(m_fields.begin(), m_fields.end(), [](const wxWeakRef<wxTextCtrl> &field) { return field.get() == nullptr; }),
                   m_fields.end());
    m_fields.emplace_back(text);
    refresh();

    auto *sizer = new wxBoxSizer(wxVERTICAL);
    sizer->Add(text, 0, wxEXPAND);
    return sizer;
}

void PreciseSeamRulesEditor::commit_pending()
{
    if (!m_edited.valid())
        return;
    // Cleared first: the snapshot and the update below may move the focus and call back here.
    const ObjectID    volume = m_edited;
    const std::string draft  = std::move(m_draft);
    m_edited                 = ObjectID();
    m_draft.clear();
    apply_precise_seam_rules(volume, draft);
}

void PreciseSeamRulesEditor::reset_from_model()
{
    m_edited = ObjectID();
    m_draft.clear();
    refresh();
}

void PreciseSeamRulesEditor::refresh()
{
    if (m_edited.valid() && m_edited == m_volume) {
        show(m_draft);
        return;
    }
    const ModelVolume *volume = m_volume.valid() ? find_volume(m_volume) : nullptr;
    show(volume ? volume->precise_seam_rules : std::string());
}

void PreciseSeamRulesEditor::rescale()
{
    for (wxTextCtrl *field : m_fields)
        if (field != nullptr)
            fit(field);
}

void PreciseSeamRulesEditor::fit(wxTextCtrl *text) const
{
    // Three lines of the current font: long rules wrap and scroll instead of growing the page.
    text->SetMinSize(wxSize(-1, 3 * text->GetCharHeight() + text->FromDIP(10)));
    wxGetApp().UpdateDarkUI(text);
}

void PreciseSeamRulesEditor::show(const std::string &text)
{
    const wxString value = from_u8(text);
    m_showing            = true;
    for (wxTextCtrl *field : m_fields)
        // An equal value keeps the cursor and the selection.
        if (field != nullptr && field->GetValue() != value)
            field->ChangeValue(value);
    m_showing = false;
}

void commit_precise_seam_rules_edit()
{
    if (PreciseSeamRulesEditor *editor = part_editor())
        editor->commit_pending();
}

void reset_precise_seam_rules_edit()
{
    if (PreciseSeamRulesEditor *editor = part_editor())
        editor->reset_from_model();
}

bool precise_seam_rules_edit_pending()
{
    const PreciseSeamRulesEditor *editor = part_editor();
    return editor != nullptr && editor->pending();
}

} // namespace GUI
} // namespace Slic3r
