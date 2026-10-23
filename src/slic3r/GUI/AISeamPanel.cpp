#include "AISeamPanel.hpp"
#include <string>
#include <functional>
#include <utility>

#include <wx/button.h>
#include <wx/dirdlg.h>
#include <wx/event.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>

#include "GUI.hpp"
#include "GUI_App.hpp"
#include "I18N.hpp"
#include "Plater.hpp"
#include "Widgets/CheckBox.hpp"
#include "libslic3r/PresetBundle.hpp"

namespace Slic3r::GUI {
namespace {
::CheckBox *create_checkbox(wxWindow *parent)
{
    auto *checkbox = new ::CheckBox(parent);
#ifndef __WXOSX__
    // Same as Field: without it MSW paints a native button face behind the checkbox bitmap.
    checkbox->SetBackgroundStyle(wxBG_STYLE_PAINT);
#endif
    return checkbox;
}

// Backend and workspace messages carry their own prefixes; the status line already names AI Seam.
wxString short_error(const std::string &error)
{
    wxString text = from_u8(error).BeforeFirst('\n');
    for (bool stripped = true; stripped;) {
        stripped = false;
        for (const wxString prefix : {wxString("AI Seam Placement: "), wxString("AI Seam: ")}) {
            wxString rest;
            if (text.StartsWith(prefix, &rest)) { text = rest; stripped = true; }
        }
    }
    return text.Left(240);
}
} // namespace

AISeamPanel::AISeamPanel(wxWindow *parent, std::function<void()> relayout)
    : wxPanel(parent), m_timer(this), m_relayout(std::move(relayout))
{
    auto *layout = new wxBoxSizer(wxVERTICAL);
    auto *enable_row = new wxBoxSizer(wxHORIZONTAL);
    enable_row->Add(new wxStaticText(this, wxID_ANY, _L("Enable AI Seam")), 1, wxALIGN_CENTER_VERTICAL);
    m_enabled = create_checkbox(this);
    m_enabled->SetToolTip(_L("Experimental AI seam exchange. Unavailable in spiral vase mode. Settings last only for this project session."));
    enable_row->Add(m_enabled, 0, wxALIGN_CENTER_VERTICAL);
    layout->Add(enable_row, 0, wxEXPAND | wxBOTTOM, 4);
    m_details = new wxPanel(this);
    auto *details = new wxBoxSizer(wxVERTICAL);
    details->Add(new wxStaticText(m_details, wxID_ANY, _L("Working folder")), 0, wxTOP, 4);
    auto *folder_row = new wxBoxSizer(wxHORIZONTAL);
    m_folder = new wxTextCtrl(m_details, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize(140, -1), wxTE_PROCESS_ENTER);
    folder_row->Add(m_folder, 1, wxALIGN_CENTER_VERTICAL);
    auto *browse = new wxButton(m_details, wxID_ANY, _L("Browse…"), wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);
    auto *open = new wxButton(m_details, wxID_ANY, _L("Open folder"), wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);
    folder_row->Add(browse, 0, wxLEFT, 4);
    folder_row->Add(open, 0, wxLEFT, 4);
    details->Add(folder_row, 0, wxEXPAND | wxBOTTOM, 6);
    auto add_toggle = [&](const wxString &label, ::CheckBox *&checkbox) {
        auto *row = new wxBoxSizer(wxHORIZONTAL);
        checkbox = create_checkbox(m_details);
        row->Add(new wxStaticText(m_details, wxID_ANY, label), 1, wxALIGN_CENTER_VERTICAL | wxRIGHT, 6);
        row->Add(checkbox, 0, wxALIGN_CENTER_VERTICAL);
        details->Add(row, 0, wxEXPAND | wxBOTTOM, 6);
        checkbox->Bind(wxEVT_TOGGLEBUTTON, [this](wxCommandEvent &event) { save(); event.Skip(); });
    };
    add_toggle(_L("Apply AI response"), m_apply);
    add_toggle(_L("Share file paths"), m_share);
    m_apply->SetToolTip(_L("Turn off to use standard seams and export a new request when slicing. The response file is kept."));
    m_share->SetToolTip(_L("Include local executable, resources, project and exchange paths in the next request. This does not grant permission to run tools."));
    // Long diagnostics wrap within the sidebar; their complete text remains in the tooltip.
    m_status = new wxStaticText(m_details, wxID_ANY, wxEmptyString);
    m_status->SetMinSize(wxSize(260, -1));
    details->Add(m_status, 0, wxEXPAND);
    m_details->SetSizer(details);
    layout->Add(m_details, 0, wxEXPAND);
    SetSizerAndFit(layout);
    m_enabled->Bind(wxEVT_TOGGLEBUTTON, [this](wxCommandEvent &event) { save(); refresh(); event.Skip(); });
    m_folder->Bind(wxEVT_TEXT_ENTER, [this](wxCommandEvent &) { save(); });
    m_folder->Bind(wxEVT_KILL_FOCUS, [this](wxFocusEvent &event) { save(); event.Skip(); });
    browse->Bind(wxEVT_BUTTON, [this](wxCommandEvent &) {
        wxDirDialog dialog(this, _L("Choose AI Seam working folder"), m_folder->GetValue(), wxDD_DIR_MUST_EXIST);
        if (dialog.ShowModal() == wxID_OK) { m_folder->SetValue(dialog.GetPath()); save(); }
    });
    open->Bind(wxEVT_BUTTON, [this](wxCommandEvent &) {
        if (!m_folder->GetValue().empty()) desktop_open_any_folder(into_u8(m_folder->GetValue()));
    });
    Bind(wxEVT_TIMER, [this](wxTimerEvent &) { refresh(); });
    wxGetApp().UpdateDarkUI(this);
    refresh();
    m_timer.Start(1000);
}

void AISeamPanel::save()
{
    if (m_syncing || !wxGetApp().plater()) return;
    AISeamPlacement::WorkspaceSettings settings;
    settings.enabled = m_enabled->GetValue();
    settings.apply_response = m_apply->GetValue();
    settings.share_paths = m_share->GetValue();
    settings.directory = into_u8(m_folder->GetValue());
    wxGetApp().plater()->set_ai_seam_settings(settings);
}

void AISeamPanel::refresh()
{
    Plater *plater = wxGetApp().plater();
    if (!plater) return;
    m_syncing = true;
    const auto settings = plater->ai_seam_settings();
    m_enabled->SetValue(settings.enabled);
    m_apply->SetValue(settings.apply_response);
    m_share->SetValue(settings.share_paths);
    if (!m_folder->HasFocus()) m_folder->ChangeValue(from_u8(settings.directory));
    const bool vase = wxGetApp().preset_bundle->prints.get_edited_preset().config.opt_bool("spiral_mode");
    m_enabled->Enable(!vase);
    const bool changed_visibility = m_details->IsShown() != settings.enabled;
    m_details->Show(settings.enabled);
    const auto result = plater->ai_seam_result();
    wxString text;
    using State = AISeamPlacement::WorkspaceState;
    switch (result.state) {
    case State::Checking: text = _L("Checking AI seam response…"); break;
    case State::Exported: text = _L("Request exported. Waiting for an AI response."); break;
    case State::Ignored: text = _L("Response ignored. New request exported."); break;
    case State::Applied:
        text = wxString::Format(_L("AI seam response applied. Seam count: %s."), from_u8(std::to_string(result.seam_count))); break;
    case State::GeometryChanged: text = _L("Geometry changed. Standard seams used; new request exported."); break;
    case State::Error:
        text = result.error.empty() ? _L("AI Seam failed. See the slicing error notification.") :
            wxString::Format(_L("AI Seam failed: %s"), short_error(result.error));
        break;
    case State::Cancelled: text = _L("AI Seam cancelled. Slice again to retry."); break;
    default: text = settings.directory.empty() ? _L("Choose a working folder.") :
        _L("Slice to export a request or apply an AI response."); break;
    }
    const bool changed_text = m_status_text != text;
    if (changed_text) {
        m_status_text = text;
        m_status->SetLabel(text);
        m_status->Wrap(260);
        m_status->SetToolTip(result.error.empty() ? text : from_u8(result.error));
        Layout();
    }
    const int plate = plater->get_partplate_list().get_curr_plate_index();
    m_folder->SetToolTip(from_u8(AISeamPlacement::workspace_filename(plate, false) + "\n" +
        AISeamPlacement::workspace_filename(plate, true)));
    if (changed_visibility || changed_text) {
        // OG_CustomCtrl measures actual child sizes, not their sizers' requested sizes.
        SetMinSize(wxSize(300, -1)); // Keep checkbox alignment stable while allowing the folded height to shrink.
        InvalidateBestSize();
        Fit();
        Layout();
        CallAfter([this] { if (m_relayout) m_relayout(); });
    }
    m_syncing = false;
}
} // namespace Slic3r::GUI
