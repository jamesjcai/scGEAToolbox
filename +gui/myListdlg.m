function [indx, tf] = myListdlg(parentfig, options, Title, ...
        prefersel, allowmulti)

    if nargin < 5, allowmulti = true; end
    if nargin < 4, prefersel = []; end

    if length(options) > 1e4
        [indx, tf] = gui.myTabledlg(parentfig, options, Title, prefersel, allowmulti);
        return; 
    end

    parentPos = parentfig.Position;
    parentCenter = [parentPos(1) + parentPos(3)/2, parentPos(2) + parentPos(4)/2];

    % Dialog size
    dlgSize = [300, 450]; % [Width, Height]

    % Compute center position
    dlgPos = [parentCenter(1) - dlgSize(1)/2, parentCenter(2) - dlgSize(2)/2, dlgSize];

    % focus(parentfig);
    % Create a modal dialog
    %    d = uifigure('Name', Title, 'Position', dlgPos, ...
    %        'WindowStyle', 'modal');

     % parentfig.WindowStyle = 'alwaysontop';
     %disp('alwaysontop')

    d = uifigure('Name', Title, 'Position', dlgPos, ...
        'WindowStyle', 'modal', 'Visible','on');
    
    if allowmulti
        multitag = 'on';
    else
        multitag = 'off';
    end

    % Create a listbox for selection
    if ~isempty(prefersel) && any(ismember(prefersel, options))
        lb = uilistbox(d, 'Items', options, 'Position', [20 60 260 370], ...
            'MultiSelect', multitag, 'Value', prefersel);
    else
        lb = uilistbox(d, 'Items', options, 'Position', [20 60 260 370], ...
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
        end
    end

     % parentfig.WindowStyle = 'normal';


    %drawnow;

    %pause(0.7);

    %d.Visible = 'on';

    % Set focus on the listbox for user interaction
    %
    % lb.focus();
    % disp('myListdlg used.');
    % Wait for user response
    %d.WindowStyle = "modal";
    uiwait(d);

    if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
        figure(parentfig);
    end

    % Get selected items
    if isvalid(d) && d.UserData
        selection = lb.Value;
        tf = 1;
        [~, indx] = ismember(selection, options);
        uiresume(d);
        delete(d);
    else
        tf = 0;
        indx = [];
        if isvalid(d)
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