function [indx, tf] = myListdlg(parentfig, options, Title, ...
    prefersel, allowmulti, allowresize, dlgSize, prompt, modal, ctxitems)
%MYLISTDLG Pick from a list, in a dialog centered on parentfig.
%
%   [indx, tf] = gui.myListdlg(parentfig, options, Title)
%   [indx, tf] = gui.myListdlg(parentfig, options, Title, prefersel, ...
%                              allowmulti, allowresize, dlgSize, prompt)
%
%   parentfig   figure to center on and return focus to; [] centers on screen
%   options     the items, a cellstr or string array
%   Title       the dialog's WINDOW title. The title bar clips anything long
%               and gives no hint that it did, so keep this to a few words
%               and put the instruction in prompt.
%   prefersel   item(s) selected initially, by value or index (default: none)
%   allowmulti  allow selecting more than one (default: true)
%   allowresize allow resizing the dialog (default: true)
%   dlgSize     [width, height] in pixels (default: [300, 450]). The height
%               is the room for the list and buttons; a prompt adds to it
%               rather than eating into the list.
%   prompt      instruction shown as a wrapped label above the list, where
%               there is room for a full sentence. '' for none, which is the
%               default and leaves the layout exactly as it was.
%   modal       block input to the parent window while the dialog is up
%               (default true). UIWAIT alone blocks only the calling code,
%               so without this the parent stays clickable and a second
%               click on the same button re-enters the callback. Pass false
%               for a dialog that must leave the parent usable.
%   ctxitems    right-click actions on the list, a struct array with fields
%               Label and Callback ([] for none, the default, which adds no
%               context menu at all). Both are called with the indices the
%               action applies to: the whole selection when the pointer is
%               over one of the selected rows, otherwise the pointed-at row
%               alone, which is also selected on the way so that the action
%               always matches what is highlighted. Right-clicking a list
%               does not move the selection by itself, which is why this is
%               done here rather than left to the caller. Label is text, or
%               a function handle of those indices returning text, so an
%               entry can name its target ("View GSE12345 on the GEO
%               website"); returning "" hides the entry. Callback gets the
%               same indices when the entry is picked. Ignored on the
%               MYTABLEDLG path below, which is taken when there are more
%               than 1e4 options.
%
%   indx is the index/indices of the chosen items, tf is 1 when OK was
%   pressed and 0 when the dialog was cancelled or closed.
%
%   Prefer Title + prompt over a long Title alone: a truncated title bar is
%   the usual reason a list dialog reads as unexplained.
%
%   See also gui.myQuestdlg, gui.myInputdlg, gui.i_centerdlgpos.

if nargin < 10, ctxitems = []; end
if nargin < 9 || isempty(modal), modal = true; end
if nargin < 8 || isempty(prompt), prompt = ''; end
if nargin < 7 || isempty(dlgSize)
    dlgSize = [300, 450]; % [Width, Height]
end
if nargin < 6, allowresize = true; end
if nargin < 5, allowmulti = true; end
if nargin < 4, prefersel = []; end

prompt = char(string(prompt));

% Reserve room for the prompt on top of the requested height, so passing one
% never shrinks the list. The line count is estimated from the dialog width
% rather than measured, which would need a realized figure; uilabel word
% wrapping does the actual breaking, so the estimate only has to be close
% enough not to clip.
promptHeight = 0;
if ~isempty(prompt)
    charsPerLine = max(20, floor((dlgSize(1)-40)/6.2));
    numLines = max(1, ceil(numel(prompt)/charsPerLine));
    promptHeight = numLines*17 + 6;
    dlgSize(2) = dlgSize(2) + promptHeight + 8;
end

% Re-entrancy guard. It still earns its place now these dialogs are modal:
% modality is switched on only once the dialog is on screen (see the note
% further down), so there is a window in which the parent is clickable,
% and a caller passing modal=false has no protection at all. Before that,
% with UIWAIT blocking only the CALLING code, clicking the same toolbar
% button built a SECOND copy of this dialog, which is what users saw on
% 'Show Cell States...' and 'Export/Save Data...'.
% GUI.I_DLGREGISTER raises the dialog that is already up; returning tf=0
% then makes callers abort on their usual cancel branch.
if ~isempty(gui.i_dlgregister(parentfig))
    indx = [];
    tf = 0;
    return;
end

% The delegation below is deliberately not registered here. GUI.MYTABLEDLG
% shares the same register and records its own dialog, so the check above
% still catches a second click while a table dialog is up.
if length(options) > 1e4
    [indx, tf] = gui.myTabledlg(parentfig, options, Title, prefersel, allowmulti);
    return;
end

% Only a window the user can see is worth focusing or raising. A hidden
% parent gets nothing: FIGURE() would show it, which is how a half-built
% GUI.MYFIGURE used to flash up empty, and FOCUS() only warns that it
% cannot focus an invisible figure.
%
% The restore on the way out goes through GUI.I_RAISEFIG, which ignores a
% handle that has since been deleted: the parent can be closed by a callback
% while the dialog is up, and a destructor that throws surfaces as a warning
% the user cannot act on.
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    if isa(parentfig, 'matlab.ui.Figure')
        try
            focus(parentfig);
            cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
        catch
            % focus() may not exist on older MATLAB; parent is brought up implicitly
        end
    else
        figure(parentfig);
        cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
    end
end


dlgPos = round(gui.i_centerdlgpos(parentfig, dlgSize));

% focus(parentfig);
% Create a modal dialog
%    d = uifigure('Name', Title, 'Position', dlgPos, ...
%        'WindowStyle', 'modal');

% parentfig.WindowStyle = 'alwaysontop';
% disp('alwaysontop')

% WindowStyle='modal' is off by default: on multi-monitor setups where the
% secondary monitor has a different DPI, MATLAB's modal centering logic uses
% an internal coordinate space that differs from MonitorPositions, causing
% the dialog to be re-centered onto the primary monitor regardless of the
% Position we set.  uiwait(d) below still blocks the calling code, so the
% dialog is functionally modal even without it.
%
% Modality is therefore switched on further down, AFTER the dialog is on
% screen, which is what GUI.MYTABLEDLG already does for the same reason:
% the re-centering happens when a figure that is not yet realized is made
% modal, so placing it first and blocking afterwards keeps both.
d = uifigure('Name', Title, 'Position', dlgPos, ...
    'Visible', 'off', 'Resize', allowresize);

gui.i_dlgregister(parentfig, d);

% pos1 = d.Position

if allowmulti
    multitag = 'on';
else
    multitag = 'off';
end

% Normalize numeric prefersel to string
if isnumeric(prefersel) && ~isempty(prefersel)
    idx = prefersel(prefersel >= 1 & prefersel <= numel(options));
    prefersel = options(idx);
end

% The prompt sits at the top; the list takes what is left between it and the
% button row. With no prompt topUsed is 20, leaving the same dlgSize(2)-60-20
% the list has always had.
topUsed = 20;
if promptHeight > 0
    uilabel(d, 'Text', prompt, 'WordWrap', 'on', ...
        'VerticalAlignment', 'top', ...
        'Position', [20, dlgSize(2)-20-promptHeight, dlgSize(1)-40, promptHeight]);
    topUsed = 20 + promptHeight + 8;
end

% Create a listbox for selection
lbHeight = dlgSize(2) - topUsed - 60;   % minus prompt, top padding, button row
if ~isempty(prefersel) && any(ismember(prefersel, options))
    lb = uilistbox(d, 'Items', options, 'Position', [20 60 dlgSize(1)-40 lbHeight], ...
        'MultiSelect', multitag, 'Value', prefersel);
else
    lb = uilistbox(d, 'Items', options, 'Position', [20 60 dlgSize(1)-40 lbHeight], ...
        'MultiSelect', multitag);
end

% Right-click actions, when the caller asked for any. The entries are made
% empty and filled on the way open, because which row the pointer is over is
% only knowable then.
if ~isempty(ctxitems)
    cm = uicontextmenu(d);
    ctxhandles = gobjects(numel(ctxitems), 1);
    for k = 1:numel(ctxitems)
        ctxhandles(k) = uimenu(cm);
    end
    cm.ContextMenuOpeningFcn = @(~, event) fillContextMenu(lb, ctxhandles, ctxitems, event);
    lb.ContextMenu = cm;
end

d.KeyPressFcn = @(src, event) jumpToFirstMatch(lb, event);

% Use UserData to track whether OK was confirmed
d.UserData = false;

% Create OK and Cancel buttons. The handles are not kept: nothing below
% refers to them, and holding them only drew a Code Analyzer warning.
uibutton(d, 'Text', 'OK', 'Position', [60 20 80 30], ...
    'ButtonPushedFcn', @(btn,event) okCallback(d));
uibutton(d, 'Text', 'Cancel', 'Position', [160 20 80 30], ...
    'ButtonPushedFcn', @(btn,event) uiresume(d));

if ~isMATLABReleaseOlderThan('R2025a')
    try
        theme(d, parentfig.Theme.BaseColorStyle);
    catch
        % theme() may not exist or parent has no Theme property; skip styling
    end
end

% d.UserData.LastState = "normal";
if ~allowresize
    d.AutoResizeChildren = 'off';
    d.SizeChangedFcn = @(src,~) enforceNormalState(src);
end

% parentfig.WindowStyle = 'normal';

% drawnow;
% pause(0.7);

% pos2 = d.Position
% assert(equal(pos1, pos2))

d.Visible = 'on';
if modal
    % Block the parent only now the dialog is placed and on screen. The
    % gap where the parent is still clickable is the one GUI.I_DLGREGISTER
    % covers, which is why that guard stays even though these are modal.
    d.WindowStyle = 'modal';
    drawnow;
    % Belt and braces against the re-centering described above: put it back
    % where i_centerdlgpos asked for. A no-op when nothing moved it.
    if ~isequal(round(d.Position), round(dlgPos))
        d.Position = dlgPos;
    end
end

% Set focus on the listbox for user interaction
%
% lb.focus();
% disp('myListdlg used.');
% Wait for user response
% d.WindowStyle = "modal";
uiwait(d);

% Get selected items
if pkg.i_isvalid(d) && d.UserData
    selection = lb.Value;
    tf = 1;
    [~, indx] = ismember(selection, options);
    uiresume(d);
    delete(d);
else
    tf = 0;
    indx = [];
    if pkg.i_isvalid(d)
        uiresume(d);
        delete(d);
    end
end

%{
Example usage:
options = {'Apple', 'Banana', 'Cherry', 'Date'};
selectedItems = gui.ui_listdlg(options, 'Select a Fruit');
disp('Selected:');
disp(selectedItems);
%}
end

function enforceNormalState(fig)
% disp('If user tries to minimize, restore immediately');

if fig.WindowState == "minimized"
    drawnow limitrate
    fig.WindowState = "normal";
end
end

function okCallback(d)
d.UserData = true;
uiresume(d);
end

function fillContextMenu(lb, handles, ctxitems, event)
% Aim every entry at the rows the action should apply to.
%   EVENT.INTERACTIONINFORMATION.ITEM is the index under the pointer, and []
%   when the click lands in the empty space below the last item - in which
%   case there is nothing to act on and every entry hides itself.
%
%   Right-clicking a uilistbox does not move the selection, so the pointer
%   and the highlight can disagree. That is resolved the way a file manager
%   does: point inside the selection and the action takes the whole of it,
%   point outside and it takes that row alone and selects it, so what runs
%   is always what is highlighted.

idx = [];
try
    info = event.InteractionInformation;
    if ~isempty(info) && isprop(info, 'Item')
        idx = info.Item;
    end
catch
    % A release that does not report what was under the pointer leaves IDX
    % empty, which hides the entries below
end

if isempty(idx)
    set(handles, 'Visible', 'off');
    return;
end

[~, sel] = ismember(string(lb.Value), string(lb.Items));
sel = sel(sel > 0);
if ismember(idx, sel)
    target = reshape(sort(sel), 1, []);   % list order, not click order
else
    target = idx;
    lb.Value = lb.Items(idx);
end

for k = 1:numel(handles)
    label = ctxitems(k).Label;
    if isa(label, 'function_handle')
        label = label(target);
    end
    label = string(label);
    if strlength(label) == 0
        handles(k).Visible = 'off';
        continue;
    end
    fcn = ctxitems(k).Callback;
    handles(k).Text = char(label);
    handles(k).MenuSelectedFcn = @(~,~) fcn(target);
    handles(k).Visible = 'on';
end
end

function jumpToFirstMatch(lb, event)
% Jump to the first item starting with the pressed letter
key = event.Character;
if isempty(key) || ~ischar(key), return; end  % Ignore non-character keys

options = lb.Items;
idx = find(startsWith(options, key, 'IgnoreCase', true), 1);
if ~isempty(idx)
    lb.Value = options{idx};  % Select matched item
end
end

%{
    fig = uifigure('Name', 'My UI App', 'Position', [500, 300, 400, 250]);

    % Create a button
    btn = uibutton(fig, 'push', ...
                   'Text', 'Click Me', ...
                   'Position', [150, 100, 100, 50], ...
                   'ButtonPushedFcn', @buttonCallback);
end

% Callback function for button press
function buttonCallback(src, event)
% uialert(src.Parent, 'Button Clicked!', 'Notification');
options = {'Apple', 'Banana', 'Cherry', 'Date'};
selectedItems = gui.ui_listdlg(options, 'Select a Fruit', src.Parent);
disp('Selected:');
disp(selectedItems);

end
%}
