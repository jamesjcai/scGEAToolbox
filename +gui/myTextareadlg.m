function [x] = myTextareadlg(parentfig, prompts, title, defaultvals, editable)


%    x = inputdlg({'Attribute Name','Attribute Values'},...
%                  'Attribute Editor', [1 80; 15 80]);

    parentPos = parentfig.Position;
    parentCenter = [parentPos(1) + parentPos(3)/2, parentPos(2) + parentPos(4)/2];

    % Dialog size
    dlgSize = [450, 300]; % [Width, Height]

    % Compute center position
    dlgPos = [parentCenter(1) - dlgSize(1)/2, parentCenter(2) - dlgSize(2)/2, dlgSize];


    % MYUIINPUTDIALOG Create an input dialog using uifigure and uitextarea
    % Similar to inputdlg but using modern UI components
    %
    % Inputs:
    %   prompts - Cell array of prompt strings
    %   title - Dialog title
    %   defaultvals - Cell array of default values (optional)
    % 
    % Output:
    %   x - Cell array of user inputs, or empty if canceled
    
    % Handle optional inputs
    if nargin < 3
        defaultvals = cell(size(prompts));
    end
    if nargin<5, editable=true(length(prompts),1); end
    
    % Create the figure
    fig = uifigure('Name', title, 'Position', ...
        dlgPos, 'Resize', 'on');
    
    % Initialize result
    x = {};
    
    % Calculate positions
    numPrompts = length(prompts);
    fieldHeight = 25;  % Height for label and single-line field
    textareaHeight = 200; % Height for multi-line textarea
    spacing = 10;      % Spacing between elements
    
    % Total height needed
    totalHeight = 0;
    for i = 1:numPrompts
        totalHeight = totalHeight + fieldHeight;
        if i == numPrompts
            totalHeight = totalHeight + textareaHeight;
        end
    end
    totalHeight = totalHeight + (numPrompts+2)*spacing;
    
    % Adjust figure height if needed
    fig.Position(4) = max(fig.Position(4), totalHeight + 60);
    
    % Width of the controls
    labelWidth = 120;
    fieldWidth = fig.Position(3) - 2*spacing - labelWidth;
    
    % Create UI components (from bottom to top)
    y = spacing;
    
    % Create OK and Cancel buttons at the bottom
    btnCancel = uibutton(fig, 'Position', [fig.Position(3)-spacing-100 y 100 30], ...
                        'Text', 'Cancel', ...
                        'ButtonPushedFcn', @(btn,event) onCancel(fig));
    
    btnOK = uibutton(fig, 'Position', [fig.Position(3)-spacing-210 y 100 30], ...
                     'Text', 'OK', ...
                     'ButtonPushedFcn', @(btn,event) onOK(fig));
    
    y = y + spacing + 30;
    
    % Create input fields (starting with the multi-line textarea at the bottom)
    inputControls = cell(numPrompts, 1);
    
    for i = numPrompts:-1:1
        % Add label
        uilabel(fig, 'Position', [spacing y+fieldHeight-20 labelWidth 20], ...
                'Text', prompts{i});
        
        % For the last field (typically "Attribute Values"), use uitextarea
        if i == numPrompts
            inputControls{i} = uitextarea(fig, ...
                'Position', [spacing+labelWidth y fieldWidth textareaHeight], ...
                'Value', defaultvals{i});
            y = y + textareaHeight + spacing;
        else
            % For other fields, use uieditfield
            inputControls{i} = uieditfield(fig, 'text', ...
                'Position', [spacing+labelWidth y fieldWidth fieldHeight], ...
                'Value', defaultvals{i});            
            y = y + fieldHeight + spacing;
        end
        if ~editable(i)
           inputControls{i}.Editable = 'off';
        end
    end

    % Make the dialog modal
    fig.WindowStyle = 'modal';
    
    % Wait for user response
    uiwait(fig);
    
    % Callback for OK button
    function onOK(fig)
        % Get values from input fields
        x = cell(numPrompts, 1);
        for j = 1:numPrompts
            x{j} = inputControls{j}.Value;
        end
        uiresume(fig);
        delete(fig);
    end
    
    % Callback for Cancel button
    function onCancel(fig)
        x = {};
        uiresume(fig);
        delete(fig);
    end
end