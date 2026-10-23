#pragma once

#include <functional>

#include <wx/panel.h>
#include <wx/timer.h>
#include <wx/string.h>

class CheckBox;
class wxStaticText;
class wxTextCtrl;

namespace Slic3r::GUI {
// Session controls share the standard options group, without adding serializable config keys.
class AISeamPanel final : public wxPanel {
public:
    AISeamPanel(wxWindow *parent, std::function<void()> relayout);
private:
    CheckBox *m_enabled = nullptr;
    CheckBox *m_apply = nullptr;
    CheckBox *m_share = nullptr;
    wxTextCtrl *m_folder = nullptr;
    wxStaticText *m_status = nullptr;
    wxPanel *m_details = nullptr;
    wxTimer m_timer;
    bool m_syncing = false;
    wxString m_status_text;
    std::function<void()> m_relayout;
    void save();
    void refresh();
};
} // namespace Slic3r::GUI
