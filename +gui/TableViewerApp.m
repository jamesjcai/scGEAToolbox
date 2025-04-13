function TableViewerApp(T, parentfig)

if nargin<2, parentfig = []; end
if nargin<1, T = []; end
    % Create a new figure window
    fig = uifigure('Name', 'Table Viewer', ...
                   'Position', [0, 0, 900, 700]);    

    try
       if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
            [px_new] = gui.i_getchildpos(parentfig, fig);
            if ~isempty(px_new)
                movegui(fig, px_new);
            else
                movegui(fig, 'center');
            end
        else
            movegui(fig, 'center');
        end
    catch
        movegui(fig, 'center');
    end

    if isempty(T)
        % Create a table with sample data
        data = generateSampleData();
    else
        T = convertvars(T, @isstring, 'cellstr');
        data = table2cell(T);
    end
    
    % Create a grid layout for better GUI organization
    mainLayout = uigridlayout(fig, [3, 1]);
    mainLayout.RowHeight = {30, '1x', 100};
    mainLayout.Padding = [10 10 10 10];
    mainLayout.RowSpacing = 10;
    
    % Create top control panel
    topPanel = uipanel(mainLayout);
    topPanel.Layout.Row = 1;
    topPanel.Layout.Column = 1;
    
    % Create top control layout
    topLayout = uigridlayout(topPanel, [1, 5]);
    topLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};
    topLayout.Padding = [5 5 5 5];
    
    % Create the table in the middle panel
    uitTable = uitable(mainLayout);
    uitTable.Data = data;
    uitTable.ColumnName = T.Properties.VariableNames;
    %{'Name', 'Age', 'Height (cm)', 'Weight (kg)', 'BMI'};
    uitTable.ColumnEditable = true;
    uitTable.Layout.Row = 2;
    uitTable.Layout.Column = 1;
    uitTable.FontSize = 12;
    
    % Enable sorting
    uitTable.ColumnSortable = true;

    % Add search functionality
    uilabel(topLayout, 'Text', 'Search:');
    searchField = uieditfield(topLayout, 'ValueChangedFcn', @(src, event) searchTable(src, uitTable));
    
    % Add row count label
    rowCountLabel = uilabel(topLayout, 'Text', sprintf('Rows: %d', size(data, 1)));
    
    % Add sort controls
    sortByDropdown = uidropdown(topLayout, 'Items', {'Select Column'}, ...
                               'ValueChangedFcn', @(src, event) sortTable(src, uitTable));
    
    % Add refresh button in top panel
    uibutton(topLayout, 'Text', 'Refresh Data', ...
             'ButtonPushedFcn', @(btn, event) refreshData(uitTable, rowCountLabel, sortByDropdown));
    

    
    % Create context menu for the table
    cm = uicontextmenu(fig);
    uitTable.ContextMenu = cm;
    
    % Add context menu items
    uimenu(cm, 'Text', 'Copy Selected Cells', 'MenuSelectedFcn', @(src, event) copySelectedCells(uitTable));
    uimenu(cm, 'Text', 'Delete Selected Rows', 'MenuSelectedFcn', @(src, event) deleteSelectedRows(uitTable, rowCountLabel));
    uimenu(cm, 'Text', 'Insert Row', 'MenuSelectedFcn', @(src, event) insertRow(uitTable, rowCountLabel));
    
    % Create a panel for the export buttons
    btnPanel = uipanel(mainLayout);
    btnPanel.Layout.Row = 3;
    btnPanel.Layout.Column = 1;
    
    % Create button layout within panel
    btnLayout = uigridlayout(btnPanel, [2, 3]);
    btnLayout.Padding = [10 10 10 10];
    btnLayout.ColumnSpacing = 20;
    btnLayout.RowSpacing = 10;
    
    % --- First Row of Buttons ---
    % Create export buttons
    exportCSVBtn = uibutton(btnLayout, 'Text', 'Export to CSV', ...
                        'ButtonPushedFcn', @(btn,event) exportToCSV(uitTable));
    exportCSVBtn.Layout.Row = 1;
    exportCSVBtn.Layout.Column = 1;
    
    exportExcelBtn = uibutton(btnLayout, 'Text', 'Export to Excel', ...
                         'ButtonPushedFcn', @(btn,event) exportToExcel(uitTable));
    exportExcelBtn.Layout.Row = 1;
    exportExcelBtn.Layout.Column = 2;
    
    exportMATBtn = uibutton(btnLayout, 'Text', 'Export to MAT File', ...
                         'ButtonPushedFcn', @(btn,event) exportToMAT(uitTable));
    exportMATBtn.Layout.Row = 1;
    exportMATBtn.Layout.Column = 3;
    
    % --- Second Row of Buttons ---
    exportWorkspaceBtn = uibutton(btnLayout, 'Text', 'Export to Workspace', ...
                              'ButtonPushedFcn', @(btn,event) exportToWorkspace(uitTable));
    exportWorkspaceBtn.Layout.Row = 2;
    exportWorkspaceBtn.Layout.Column = 1;
    
    printBtn = uibutton(btnLayout, 'Text', 'Print Table', ...
                   'ButtonPushedFcn', @(btn,event) printTable(uitTable));
    printBtn.Layout.Row = 2;
    printBtn.Layout.Column = 2;
    
    statsBtn = uibutton(btnLayout, 'Text', 'Show Statistics', ...
                    'ButtonPushedFcn', @(btn,event) showStatistics(uitTable));
    statsBtn.Layout.Row = 2;
    statsBtn.Layout.Column = 3;
    
    
    % Update sort dropdown with column names
    updateSortDropdown(sortByDropdown, uitTable);
    
    % Create callback to update row count when table changes
    addlistener(uitTable, 'Data', 'PostSet', @(src, event) updateRowCount(rowCountLabel, uitTable));
end

% function data = generateSampleData()
%     % Generate some sample data for the table
%     names = {'John Smith', 'Mary Johnson', 'Robert Williams', 'Susan Brown', ...
%              'Michael Davis', 'Patricia Miller', 'James Wilson', 'Linda Moore', ...
%              'David Taylor', 'Jennifer Anderson', 'William Thomas', 'Elizabeth Jackson'};
%     ages = randi([18, 65], length(names), 1);
%     heights = randi([150, 190], length(names), 1);
%     weights = randi([50, 100], length(names), 1);
% 
%     % Calculate BMI (weight in kg / (height in m)^2)
%     bmi = weights ./ ((heights/100).^2);
%     bmi = round(bmi * 10) / 10;  % Round to 1 decimal place
% 
%     % Combine all data into a table
%     data = table(names', ages, heights, weights, bmi, ...
%                 'VariableNames', {'Name', 'Age', 'Height', 'Weight', 'BMI'});
%     data = table2cell(data);
% end

function exportToCSV(tableObj)
    % Function to export table data to a CSV file
    [file, path] = uiputfile('*.csv', 'Save Table Data');
    if isequal(file, 0) || isequal(path, 0)
        % User canceled the dialog
        return;
    end
    
    % Get the full file path
    fullFilePath = fullfile(path, file);
    
    % Get the column names and data
    columnNames = tableObj.ColumnName;
    data = tableObj.Data;
    
    % Convert to MATLAB table
    T = cell2table(data, 'VariableNames', strrep(columnNames, ' ', '_'));
    
    % Write to CSV
    writetable(T, fullFilePath);
    
    % Display confirmation message
    uialert(tableObj.Parent.Parent, ['Table data saved to: ' fullFilePath], 'Export Successful', 'Icon', 'success');
end

function exportToExcel(tableObj)
    % Function to export table data to an Excel file
    [file, path] = uiputfile({'*.xlsx;*.xls', 'Excel Files (*.xlsx, *.xls)'}, 'Save Table Data');
    if isequal(file, 0) || isequal(path, 0)
        % User canceled the dialog
        return;
    end
    
    % Get the full file path
    fullFilePath = fullfile(path, file);
    
    % Get the column names and data
    columnNames = tableObj.ColumnName;
    data = tableObj.Data;
    
    % Convert to MATLAB table
    T = cell2table(data, 'VariableNames', strrep(columnNames, ' ', '_'));
    
    % Write to Excel
    writetable(T, fullFilePath, 'Sheet', 'TableData', 'WriteVariableNames', true);
    
    % Display confirmation message
    uialert(tableObj.Parent.Parent, ['Table data saved to: ' fullFilePath], 'Export Successful', 'Icon', 'success');
end

function exportToMAT(tableObj)
    % Function to export table data to a MAT file
    [file, path] = uiputfile('*.mat', 'Save Table Data');
    if isequal(file, 0) || isequal(path, 0)
        % User canceled the dialog
        return;
    end
    
    % Get the full file path
    fullFilePath = fullfile(path, file);
    
    % Get the column names and data
    columnNames = tableObj.ColumnName;
    data = tableObj.Data;
    
    % Convert to MATLAB table
    tableData = cell2table(data, 'VariableNames', strrep(columnNames, ' ', '_'));
    
    % Save to MAT file
    save(fullFilePath, 'tableData');
    
    % Display confirmation message
    uialert(tableObj.Parent.Parent, ['Table data saved to: ' fullFilePath], 'Export Successful', 'Icon', 'success');
end

function exportToWorkspace(tableObj)
    % Function to export table data to the workspace
    
    % Create input dialog to get variable name
    prompt = {'Enter variable name:'};
    dlgtitle = 'Export to Workspace';
    dims = [1 40];
    definput = {'tableData'};
    
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    
    if isempty(answer)
        % User canceled the dialog
        return;
    end
    
    % Get variable name
    varName = answer{1};
    
    % Get the column names and data
    columnNames = tableObj.ColumnName;
    data = tableObj.Data;
    
    % Convert to MATLAB table
    tableData = cell2table(data, 'VariableNames', strrep(columnNames, ' ', '_'));
    
    % Assign to workspace variable
    assignin('base', varName, tableData);
    
    % Display confirmation message
    uialert(tableObj.Parent.Parent, ['Table data exported to workspace variable: ' varName], 'Export Successful', 'Icon', 'success');
end

function refreshData(tableObj, rowCountLabel, sortByDropdown)
    % Function to refresh data with new random values
    tableObj.Data = generateSampleData();
    
    % Update row count
    updateRowCount(rowCountLabel, tableObj);
    
    % Reset sort dropdown
    updateSortDropdown(sortByDropdown, tableObj);
end

function printTable(tableObj)
    % Function to print or preview the table
    
    % Create a new figure for printing
    printFig = figure('Visible', 'off', 'Name', 'Print Preview');
    
    % Get data and column names
    data = tableObj.Data;
    colNames = tableObj.ColumnName;
    
    % Create a uitable in the print figure
    printUIT = uitable(printFig);
    printUIT.Data = data;
    printUIT.ColumnName = colNames;
    
    % Position the table (this is not in a grid layout, so it's OK)
    printUIT.Units = 'normalized';
    printUIT.Position = [0.02 0.02 0.96 0.96];
    
    % Set figure title
    title = annotation(printFig, 'textbox', [0.3, 0.95, 0.4, 0.05], 'String', 'Table Data');
    title.HorizontalAlignment = 'center';
    title.FontSize = 14;
    title.FontWeight = 'bold';
    title.EdgeColor = 'none';
    
    % Show print preview
    printFig.Visible = 'on';
    printpreview(printFig);
end

function showStatistics(tableObj)
    % Function to show basic statistics for numeric columns
    
    % Get data and column names
    data = tableObj.Data;
    colNames = tableObj.ColumnName;
    
    % Convert cell array to numeric where possible
    numericData = [];
    numericColNames = {};
    
    for col = 1:size(data, 2)
        % Try to convert this column to numeric
        colData = data(:, col);
        if all(cellfun(@(x) isnumeric(x) || (ischar(x) && ~isempty(str2double(x)) && ~isnan(str2double(x))), colData))
            numericData = [numericData, cellfun(@double, colData)];
            numericColNames{end+1} = colNames{col};
        end
    end
    
    % If no numeric data, show message and return
    if isempty(numericData)
        uialert(tableObj.Parent.Parent, 'No numeric data found in the table.', 'Statistics', 'Icon', 'info');
        return;
    end
    
    % Calculate statistics
    meanVals = mean(numericData, 1);
    stdVals = std(numericData, 1);
    minVals = min(numericData, [], 1);
    maxVals = max(numericData, [], 1);
    medianVals = median(numericData, 1);
    
    % Create statistics table
    statNames = {'Mean', 'Std Dev', 'Minimum', 'Maximum', 'Median'};
    statData = [meanVals; stdVals; minVals; maxVals; medianVals];
    
    % Create statistics figure
    statsFig = uifigure('Name', 'Table Statistics', 'Position', [200, 200, 600, 400]);
    
    % Create a grid layout for the statistics
    statsLayout = uigridlayout(statsFig, [1, 1]);
    
    % Add table to figure using the grid layout
    statsTable = uitable(statsLayout);
    statsTable.Data = statData;
    statsTable.RowName = statNames;
    statsTable.ColumnName = numericColNames;
    
    % No need to set position when using grid layout
    % Adjust column widths
    statsTable.ColumnWidth = {100};
end

function searchTable(searchField, tableObj)
    % Function to search the table
    searchText = lower(searchField.Value);
    
    % If search text is empty, show all rows
    if isempty(searchText)
        tableObj.Data = generateSampleData();
        return;
    end
    
    % Get data and search
    data = tableObj.Data;
    matchRows = false(size(data, 1), 1);
    
    % Search each cell
    for row = 1:size(data, 1)
        for col = 1:size(data, 2)
            cellValue = data{row, col};
            % Convert to string if numeric
            if isnumeric(cellValue)
                cellValue = num2str(cellValue);
            elseif ~ischar(cellValue) && ~isstring(cellValue)
                continue;
            end
            
            % Check for match
            if contains(lower(cellValue), searchText)
                matchRows(row) = true;
                break;
            end
        end
    end
    
    % Update table with matching rows
    tableObj.Data = data(matchRows, :);
end

function updateSortDropdown(dropdown, tableObj)
    % Function to update sort dropdown items with column names
    % Convert column names to cell array if it's not already
    if ~iscell(tableObj.ColumnName)
        columnNames = cellstr(tableObj.ColumnName);
    else
        columnNames = tableObj.ColumnName;
    end
    
    % Ensure Select Column is a cell string
    selectColumn = {'Select Column'};
    
    % Combine items ensuring they are compatible types
    items = [selectColumn; columnNames];
    
    dropdown.Items = items;
    dropdown.Value = 'Select Column';
end

function sortTable(dropdown, tableObj)
    % Function to sort the table by selected column
    selectedColumn = dropdown.Value;
    
    % If "Select Column" is selected, do nothing
    if strcmp(selectedColumn, 'Select Column')
        return;
    end
    
    % Get data and column index
    data = tableObj.Data;
    colIndex = strcmp(tableObj.ColumnName, selectedColumn);
    
    % Extract the column to sort
    colData = data(:, colIndex);
    
    % Convert all cells to string for sorting
    strColData = cell(size(colData));
    for i = 1:length(colData)
        if isnumeric(colData{i})
            strColData{i} = num2str(colData{i});
        else
            strColData{i} = colData{i};
        end
    end
    
    % Sort
    [~, sortIdx] = sort(strColData);
    
    % Update table
    tableObj.Data = data(sortIdx, :);
end

function updateRowCount(label, tableObj)
    % Function to update row count label
    label.Text = sprintf('Rows: %d', size(tableObj.Data, 1));
end

function copySelectedCells(tableObj)
    % Function to copy selected cells to clipboard
    if isempty(tableObj.Selection)
        return;
    end
    
    % Get selected cells
    selectedRows = unique(tableObj.Selection(:, 1));
    selectedCols = unique(tableObj.Selection(:, 2));
    
    % Extract data
    data = tableObj.Data;
    selectedData = data(selectedRows, selectedCols);
    
    % Create string for clipboard
    clipboardStr = '';
    for row = 1:size(selectedData, 1)
        rowStr = '';
        for col = 1:size(selectedData, 2)
            cellValue = selectedData{row, col};
            if isnumeric(cellValue)
                cellValue = num2str(cellValue);
            end
            rowStr = [rowStr, cellValue, '\t'];
        end
        clipboardStr = [clipboardStr, rowStr(1:end-2), '\n'];
    end
    
    % Copy to clipboard
    clipboard('copy', clipboardStr);
end

function deleteSelectedRows(tableObj, rowCountLabel)
    % Function to delete selected rows
    if isempty(tableObj.Selection)
        return;
    end
    
    % Get selected rows
    selectedRows = unique(tableObj.Selection(:, 1));
    
    % Get data
    data = tableObj.Data;
    
    % Create mask of rows to keep
    keepRows = true(size(data, 1), 1);
    keepRows(selectedRows) = false;
    
    % Update table
    tableObj.Data = data(keepRows, :);
    
    % Update row count
    updateRowCount(rowCountLabel, tableObj);
end

function insertRow(tableObj, rowCountLabel)
    % Function to insert a new row
    data = tableObj.Data;
    
    % Create empty row
    emptyRow = cell(1, size(data, 2));
    
    % Get selected row or add at end
    if isempty(tableObj.Selection)
        % Add at end
        data = [data; emptyRow];
    else
        % Insert after selected row
        selectedRow = tableObj.Selection(1, 1);
        data = [data(1:selectedRow, :); emptyRow; data(selectedRow+1:end, :)];
    end
    
    % Update table
    tableObj.Data = data;
    
    % Update row count
    updateRowCount(rowCountLabel, tableObj);
end