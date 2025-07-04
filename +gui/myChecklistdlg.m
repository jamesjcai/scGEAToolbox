function selectedItems = myChecklistdlg(parentfig, items, varargin)
    % CHECKLISTDIALOG Create a dialog with checklist for multi-selection
    %
    % selectedItems = checklistDialog(items) creates a dialog with a checklist
    % of items for user selection. Returns cell array of selected items.
    %
    % selectedItems = checklistDialog(items, 'Title', title) sets dialog title
    % selectedItems = checklistDialog(items, 'Prompt', prompt) sets prompt text
    % selectedItems = checklistDialog(items, 'DefaultSelection', indices) 
    %   sets default selected items by indices
    %
    % Example:
    %   items = {'Option 1', 'Option 2', 'Option 3', 'Option 4'};
    %   selected = gui.myChecklistdlg(items, 'Title', 'Select Options');
    
    % Parse input arguments
    p = inputParser;
    addParameter(p, 'parentfig', []);
    addRequired(p, 'items', @(x) iscell(x) || isstring(x) || ischar(x));
    addParameter(p, 'Title', 'Select Items', @ischar);
    addParameter(p, 'Prompt', 'Select one or more items:', @ischar);
    addParameter(p, 'DefaultSelection', [], @isnumeric);
    parse(p, items, varargin{:});
    
    % Convert items to cell array if needed
    if ischar(items)
        items = {items};
    elseif isstring(items)
        items = cellstr(items);
    end
    
    % Initialize return variable
    selectedItems = {};
    

    dialogWidth = 350;
    dialogHeight = 400;
    dialogX = 100;
    dialogY = 100;
    if ~isempty(parentfig)
        parentPos = parentfig.Position;
        dialogX = parentPos(1) + (parentPos(3) - dialogWidth) / 2;
        dialogY = parentPos(2) + (parentPos(4) - dialogHeight) / 2;
    end

    % Create the dialog figure
    fig = uifigure('Name', p.Results.Title, ...
                   'Position', [dialogX, dialogY, dialogWidth, dialogHeight], ...
                   'WindowStyle', 'modal', ...
                   'Resize', 'off');
    
    % Create main grid layout
    mainGrid = uigridlayout(fig, [4, 1]);
    mainGrid.RowHeight = {'fit', '1x', 'fit', 'fit'};
    mainGrid.Padding = [20, 20, 20, 20];
    mainGrid.RowSpacing = 15;
    
    % Add prompt label
    promptLabel = uilabel(mainGrid, 'Text', p.Results.Prompt, ...
                         'FontSize', 12, 'FontWeight', 'bold');
    
    % Create scrollable panel for checklist
    scrollPanel = uipanel(mainGrid, 'BorderType', 'line', ...
                         'BackgroundColor', [1, 1, 1]);
    
    % Create grid layout for checkboxes
    numItems = length(items);
    checkGrid = uigridlayout(scrollPanel, [numItems, 1]);
    checkGrid.RowHeight = repmat({'fit'}, 1, numItems);
    checkGrid.Padding = [10, 10, 10, 10];
    checkGrid.RowSpacing = 5;
    
    % Create checkboxes
    checkboxes = cell(numItems, 1);
    for i = 1:numItems
        checkboxes{i} = uicheckbox(checkGrid, 'Text', items{i}, ...
                                  'FontSize', 10);
        % Set default selection if specified
        if ismember(i, p.Results.DefaultSelection)
            checkboxes{i}.Value = true;
        end
    end
    
    % Create selection info label
    selectionInfo = uilabel(mainGrid, 'Text', 'Selected: 0 items', ...
                           'FontSize', 10, 'FontColor', [0.5, 0.5, 0.5]);
    
    % Add callback to update selection count
    updateSelectionCount();
    for i = 1:numItems
        checkboxes{i}.ValueChangedFcn = @(~,~) updateSelectionCount();
    end
    
    % Create button panel
    buttonPanel = uipanel(mainGrid, 'BorderType', 'none', ...
                         'BackgroundColor', fig.Color);
    buttonGrid = uigridlayout(buttonPanel, [1, 4]);
    buttonGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    buttonGrid.ColumnSpacing = 10;
    
    % Spacer
    uilabel(buttonGrid, 'Text', '');
    
    % Select All button
    selectAllBtn = uibutton(buttonGrid, 'Text', 'Select All', ...
                           'ButtonPushedFcn', @selectAllCallback);
    
    % Clear All button
    clearAllBtn = uibutton(buttonGrid, 'Text', 'Clear All', ...
                          'ButtonPushedFcn', @clearAllCallback);
    
    % OK button
    okBtn = uibutton(buttonGrid, 'Text', 'OK', ...
                    'ButtonPushedFcn', @okCallback, ...
                    'FontWeight', 'bold');
    
    % Cancel button
    cancelBtn = uibutton(buttonGrid, 'Text', 'Cancel', ...
                        'ButtonPushedFcn', @cancelCallback);
    
    % Wait for user interaction
    uiwait(fig);
    
    % Nested callback functions
    function updateSelectionCount()
        selectedCount = sum(cellfun(@(x) x.Value, checkboxes));
        selectionInfo.Text = sprintf('Selected: %d items', selectedCount);
    end
    
    function selectAllCallback(~, ~)
        for j = 1:numItems
            checkboxes{j}.Value = true;
        end
        updateSelectionCount();
    end
    
    function clearAllCallback(~, ~)
        for j = 1:numItems
            checkboxes{j}.Value = false;
        end
        updateSelectionCount();
    end
    
    function okCallback(~, ~)
        % Get selected items
        selectedIndices = cellfun(@(x) x.Value, checkboxes);
        selectedItems = items(selectedIndices);
        
        % Close dialog
        if isvalid(fig)
            delete(fig);
        end
    end
    
    function cancelCallback(~, ~)
        % Return empty selection
        selectedItems = {};
        
        % Close dialog
        if isvalid(fig)
            delete(fig);
        end
    end
end