function [answer] = myInputdlg(prompt, dlgtitle, definput, parentfig)

    % if nargin<1, parentfig = uifigure('Visible', 'off'); end
    % 
    % figure(parentfig);
    % % Example usage with multiple inputs
    % prompt = {'Enter matrix size:', 'Enter colormap name:'};
    % dlgtitle = 'Input';
    % % fieldsize = [1 45; 1 45]; % Not used in uifigure, but kept for consistency
    % definput = {'20', 'hsv'};

    answer = openInputDialog(prompt, dlgtitle, definput, parentfig);
    
    % Display result
    if ~isempty(answer)
        disp('User Input:');
        disp(answer);
    else
        disp('User cancelled the input dialog.');
    end
end

function answer = openInputDialog(prompt, dlgtitle, definput, parentfig)

     dialogWidth = 350;
     dialogHeight = 100 + 50 * numel(prompt); % Adjust height based on the number of inputs

    % Compute center position relative to parentFig
    parentPos = parentfig.Position;
    dialogX = parentPos(1) + (parentPos(3) - dialogWidth) / 2;
    dialogY = parentPos(2) + (parentPos(4) - dialogHeight) / 2;

    % % Parent figure (invisible, to keep the UI elements modal)
    % parentFig = uifigure('Visible', 'off');
    % 
    % % Dialog size

    % 
    % % Compute center position on screen
    % screenSize = get(groot, 'ScreenSize');
    % dialogX = (screenSize(3) - dialogWidth) / 2;
    % dialogY = (screenSize(4) - dialogHeight) / 2;

    % Create modal dialog
    d = uifigure('Position', [dialogX, dialogY, dialogWidth, dialogHeight], ...
                 'Name', dlgtitle, 'WindowStyle', 'modal');

    % Store user input
    numFields = numel(prompt);
    fields = gobjects(numFields, 1);
    
    % Create UI elements
    for i = 1:numFields
        uilabel(d, 'Text', prompt{i}, ...
                   'Position', [20, dialogHeight - 40 - (i - 1) * 50, 310, 22]);
        fields(i) = uieditfield(d, 'text', ...
                                'Position', [20, dialogHeight - 65 - (i - 1) * 50, 310, 22], ...
                                'Value', definput{i});
    end

    % OK button
    okBtn = uibutton(d, 'Text', 'OK', ...
                     'Position', [100, 20, 60, 30], ...
                     'ButtonPushedFcn', @(btn, event) onOKButton(d, fields));

    % Cancel button
    cancelBtn = uibutton(d, 'Text', 'Cancel', ...
                         'Position', [190, 20, 60, 30], ...
                         'ButtonPushedFcn', @(btn, event) delete(d));

    % Wait for user input
    uiwait(d);
    
    % Retrieve data
    if isvalid(d) % If the dialog was not closed by user
        answer = arrayfun(@(f) f.Value, fields, 'UniformOutput', false);
        delete(d);
    else
        answer = {};
    end
end

function onOKButton(d, ~)
    % Resume UI execution when OK is pressed
    uiresume(d);
end
