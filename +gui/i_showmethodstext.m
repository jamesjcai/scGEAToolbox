function fig = i_showmethodstext(parentfig, txt)
%I_SHOWMETHODSTEXT Show Methods text in a window, with copy and save buttons.
%
%   gui.i_showmethodstext(parentfig, txt)
%   fig = gui.i_showmethodstext(parentfig, txt)
%
%   A read-only, wrapped view of TXT with three buttons: Copy to Clipboard,
%   Save as Text File... and Close. Nothing is copied until the user clicks
%   Copy; the status line then says it was. Modal, and centred on
%   PARENTFIG, like the toolbox's other dialogs.
%
%   With an output argument the window is returned without waiting for it
%   to close, so a test can press its buttons.
%
%   See also GUI.CALLBACK_SHOWMETHODSTEXT, SC_METHODSTEXT.

arguments
    parentfig
    txt (1,1) string
end

dlgSize = [720, 520];
fig = uifigure('Name', 'Methods Text', ...
    'Position', gui.i_centerdlgpos(parentfig, dlgSize), 'Resize', 'on');

g = uigridlayout(fig, [3, 1]);
g.RowHeight = {'fit', '1x', 'fit'};
g.Padding = [12 12 12 12];
g.RowSpacing = 8;

uilabel(g, 'WordWrap', 'on', 'Text', ['Methods text for the analysis ' ...
    'steps run on this dataset, with the settings used. Read it before ' ...
    'pasting: it says what the toolbox did, not why.']);

ta = uitextarea(g, 'Value', splitlines(txt), 'Editable', 'off', ...
    'WordWrap', 'on', 'Tag', 'MethodsText');

bar = uigridlayout(g, [1, 4]);
bar.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
bar.Padding = [0 0 0 0];
status = uilabel(bar, 'Text', '', 'Tag', 'Status');
uibutton(bar, 'Text', 'Copy to Clipboard', 'Tag', 'Copy', ...
    'ButtonPushedFcn', @(~, ~) i_copy());
uibutton(bar, 'Text', 'Save as Text File...', 'Tag', 'Save', ...
    'ButtonPushedFcn', @(~, ~) i_save());
uibutton(bar, 'Text', 'Close', 'Tag', 'Close', ...
    'ButtonPushedFcn', @(~, ~) delete(fig));

if nargout > 0
    return
end
fig.WindowStyle = 'modal';
uiwait(fig);
if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
    figure(parentfig);
end

    function i_copy()
        clipboard('copy', strjoin(string(ta.Value), newline));
        status.Text = 'Copied to the clipboard.';
    end

    function i_save()
        [file, path] = uiputfile('*.txt', 'Save Methods Text', 'methods_text.txt');
        figure(fig);
        if isequal(file, 0), return; end
        fid = fopen(fullfile(path, file), 'w', 'n', 'UTF-8');
        if fid < 0
            status.Text = 'Could not write that file.';
            return
        end
        closer = onCleanup(@() fclose(fid));
        fprintf(fid, '%s\n', strjoin(string(ta.Value), newline));
        status.Text = sprintf('Saved to %s.', file);
    end
end
