function callback_ShortcutsGuide(src, ~)
%CALLBACK_SHORTCUTSGUIDE List the keyboard shortcuts of the app's menus.
%
%   Behind Help > Keyboard Shortcuts. The list is read from the menus'
%   own Accelerator properties rather than from a written table: the table
%   it replaces (a row of assets/Misc/refinfo.txt) had gone stale, naming
%   items under old labels and missing Ctrl+Z, and a list built from the
%   menus cannot drift from them.
%
%   See also GUI.GUI_SHOWREFINFO.

[FigureHandle, ~] = gui.gui_getfigsce(src);
lines = gui.i_menushortcuts(FigureHandle);
if isempty(lines)
    gui.myHelpdlg(FigureHandle, 'No menu item has a keyboard shortcut.', 'Shortcuts');
    return;
end
gui.myHelpdlg(FigureHandle, strjoin(lines, newline), 'Keyboard Shortcuts');
end
