function [indx, tf] = myListdlg(parentFig, options, Title)

    parentPos = parentFig.Position;
    parentCenter = [parentPos(1) + parentPos(3)/2, parentPos(2) + parentPos(4)/2];

    % Dialog size
    dlgSize = [300, 450]; % [Width, Height]

    % Compute center position
    dlgPos = [parentCenter(1) - dlgSize(1)/2, parentCenter(2) - dlgSize(2)/2, dlgSize];

    % Create a modal dialog
    d = uifigure('Name', Title, 'Position', dlgPos, ...
        'WindowStyle', 'modal');


    % Create a listbox for selection
    lb = uilistbox(d, 'Items', options, 'Position', [20 60 260 370], 'MultiSelect', 'on');

    % Create OK button
    btnOK = uibutton(d, 'Text', 'OK', 'Position', [60 20 80 30], ...
        'ButtonPushedFcn', @(btn,event) uiresume(d));

    % Create Cancel button
    btnCancel = uibutton(d, 'Text', 'Cancel', 'Position', [160 20 80 30], ...
        'ButtonPushedFcn', @(btn,event) (close(d)));
    % Wait for user response
    uiwait(d);
    
    % Get selected items
    if isvalid(lb)
        selection = lb.Value;
        delete(d);
        tf = 1;
        [~, indx]=ismember(selection, options);
    else
        selection = {};
        tf = 0;
        indx = [];
    end

    %{
    Example usage:
    options = {'Apple', 'Banana', 'Cherry', 'Date'};
    selectedItems = gui.ui_listdlg(options, 'Select a Fruit');
    disp('Selected:');
    disp(selectedItems);
    %}
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