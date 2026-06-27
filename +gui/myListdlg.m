function [indx, tf] = myListdlg(parentfig, options, Title, ...
    prefersel, allowmulti, allowresize, dlgSize)

if nargin < 7
    dlgSize = [300, 450]; % [Width, Height]
end
if nargin < 6, allowresize = true; end
if nargin < 5, allowmulti = true; end
if nargin < 4, prefersel = []; end

if length(options) > 1e4
    [indx, tf] = gui.myTabledlg(parentfig, options, Title, prefersel, allowmulti);
    return;
end

if ~isempty(parentfig)
    if isa(parentfig, 'matlab.ui.Figure')
        try
            focus(parentfig);
            cleanupObj = onCleanup(@() focus(parentfig));
        catch
            % focus() may not exist on older MATLAB; parent is brought up implicitly
        end
    else
        figure(parentfig);
        cleanupObj = onCleanup(@() figure(parentfig));
    end
end


dlgPos = round(gui.i_centerdlgpos(parentfig, dlgSize));

% focus(parentfig);
% Create a modal dialog
%    d = uifigure('Name', Title, 'Position', dlgPos, ...
%        'WindowStyle', 'modal');

 % parentfig.WindowStyle = 'alwaysontop';
 % disp('alwaysontop')

% WindowStyle='modal' is intentionally omitted: on multi-monitor setups
% where the secondary monitor has a different DPI, MATLAB's modal centering
% logic uses an internal coordinate space that differs from MonitorPositions,
% causing the dialog to be re-centered onto the primary monitor regardless
% of the Position we set.  uiwait(d) below still blocks the calling code,
% so the dialog is functionally modal.
d = uifigure('Name', Title, 'Position', dlgPos, ...
    'Visible', 'off', 'Resize', allowresize);

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

% Create a listbox for selection
lbHeight = dlgSize(2) - 60 - 20;   % dialog height minus button row and top padding
if ~isempty(prefersel) && any(ismember(prefersel, options))
    lb = uilistbox(d, 'Items', options, 'Position', [20 60 dlgSize(1)-40 lbHeight], ...
        'MultiSelect', multitag, 'Value', prefersel);
else
    lb = uilistbox(d, 'Items', options, 'Position', [20 60 dlgSize(1)-40 lbHeight], ...
        'MultiSelect', multitag);
end

d.KeyPressFcn = @(src, event) jumpToFirstMatch(lb, event);

% Use UserData to track whether OK was confirmed
d.UserData = false;

% Create OK button
btnOK = uibutton(d, 'Text', 'OK', 'Position', [60 20 80 30], ...
'ButtonPushedFcn', @(btn,event) okCallback(d));

% Create Cancel button
btnCancel = uibutton(d, 'Text', 'Cancel', 'Position', [160 20 80 30], ...
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
