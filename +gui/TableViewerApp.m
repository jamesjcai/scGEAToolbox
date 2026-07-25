function [fig] = TableViewerApp(T, parentfig, defname, customActions)
% TABLEVIEWERAPP Display one or more tables in an interactive viewer.
%
% fig = TABLEVIEWERAPP(T, parentfig, defname, customActions) shows table T.
%
% T may also be a struct array with fields Name and Table, in which case a
% "View" dropdown lets the user switch between the named tables. For
% example, DE results can be shown as:
%   views(1) = struct('Name', 'All genes', 'Table', T);
%   views(2) = struct('Name', 'Up-regulated', 'Table', Tup);
%   views(3) = struct('Name', 'Down-regulated', 'Table', Tdn);
% Search, sorting, and the export buttons act on the displayed view.

if nargin < 4 || isempty(customActions), customActions = struct([]); end
if nargin<3, defname = []; end
if nargin<2, parentfig = []; end
if nargin<1, T = []; end

views = i_normalizeviews(T);

    % Create a new figure window
fig = uifigure('Name', 'Table Viewer', ...
                   'Position', [0, 0, 900, 480],'Visible','off');

gui.i_movegui2parent(fig, parentfig);

fig.UserData = struct('Views', views, 'CurrentIndex', 1);
[data, headers] = i_viewcontent(views(1));

    % Create a grid layout for better GUI organization
mainLayout = uigridlayout(fig, [3, 1]);
mainLayout.RowHeight = {30, '1x', 100};
mainLayout.Padding = [10 10 10 10];
mainLayout.RowSpacing = 10;

    % Create top control panel
topPanel = uipanel(mainLayout);
topPanel.Layout.Row = 1;
topPanel.Layout.Column = 1;

    % Create top control layout; a "View" selector is added only when the
    % caller supplied more than one named table.
hasViewSelector = numel(views) > 1;
nTopCols = 5 + 2*hasViewSelector;
topLayout = uigridlayout(topPanel, [1, nTopCols]);
topLayout.ColumnWidth = repmat({'1x'}, 1, nTopCols);
topLayout.Padding = [5 5 5 5];

    % Create the table in the middle panel
uitTable = uitable(mainLayout);
uitTable.Data = data;
uitTable.ColumnName = headers;
uitTable.ColumnEditable = true;
uitTable.Layout.Row = 2;
uitTable.Layout.Column = 1;
uitTable.FontSize = 12;

    % Enable sorting
uitTable.ColumnSortable = true;

    % Add view selector. Its callback needs controls created below, so it is
    % attached once they exist.
if hasViewSelector
    uilabel(topLayout, 'Text', 'View:');
    viewDropdown = uidropdown(topLayout, 'Items', {views.Name});
end

    % Add search functionality
uilabel(topLayout, 'Text', 'Search:');
searchField = uieditfield(topLayout, 'ValueChangedFcn', @(src, event) searchTable(src, uitTable, fig));

    % Add row count label
rowCountLabel = uilabel(topLayout, 'Text', sprintf('Rows: %d', size(data, 1)));

    % Add sort controls
sortByDropdown = uidropdown(topLayout, 'Items', {'Select Column'}, ...
                               'ValueChangedFcn', @(src, event) sortTable(src, uitTable));

    % Add refresh button in top panel
uibutton(topLayout, 'Text', 'Refresh Data', ...
             'ButtonPushedFcn', @(btn, event) refreshData(uitTable, rowCountLabel, sortByDropdown, fig));

if hasViewSelector
    viewDropdown.ValueChangedFcn = @(src, event) switchView(src, fig, ...
        uitTable, searchField, rowCountLabel, sortByDropdown);
end


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

secondRowCount = 1 + numel(customActions);
nBtnCols = max(3, secondRowCount);

    % Create button layout within panel
btnLayout = uigridlayout(btnPanel, [2, nBtnCols]);
btnLayout.Padding = [10 10 10 10];
btnLayout.ColumnSpacing = 20;
btnLayout.RowSpacing = 10;
btnLayout.ColumnWidth = repmat({'1x'}, 1, nBtnCols);

    % --- First Row of Buttons ---
    % Create export buttons
exportCSVBtn = uibutton(btnLayout, 'Text', 'Export to CSV', ...
                        'ButtonPushedFcn', @(btn,event) exportToCSV(uitTable, i_defname(fig, defname)));
exportCSVBtn.Layout.Row = 1;
exportCSVBtn.Layout.Column = 1;

exportExcelBtn = uibutton(btnLayout, 'Text', 'Export to Excel', ...
                         'ButtonPushedFcn', @(btn,event) exportToExcel(uitTable, i_defname(fig, defname)));
exportExcelBtn.Layout.Row = 1;
exportExcelBtn.Layout.Column = 2;

exportMATBtn = uibutton(btnLayout, 'Text', 'Export to MAT File', ...
                         'ButtonPushedFcn', @(btn,event) exportToMAT(uitTable, i_defname(fig, defname)));
exportMATBtn.Layout.Row = 1;
exportMATBtn.Layout.Column = 3;

    % --- Second Row of Buttons ---
exportWorkspaceBtn = uibutton(btnLayout, 'Text', 'Export to Workspace', ...
                              'ButtonPushedFcn', @(btn,event) exportToWorkspace(uitTable, fig));
exportWorkspaceBtn.Layout.Row = 2;
exportWorkspaceBtn.Layout.Column = 1;

for actionIdx = 1:numel(customActions)
    action = customActions(actionIdx);
    actionBtn = uibutton(btnLayout, 'Text', action.Text, ...
        'ButtonPushedFcn', @(btn, event) executeCustomAction(action, uitTable, fig));
    actionBtn.Layout.Row = 2;
    actionBtn.Layout.Column = actionIdx + 1;
    if isfield(action, 'Tooltip') && ~isempty(action.Tooltip)
        actionBtn.Tooltip = action.Tooltip;
    end
end

    %{
    printBtn = uibutton(btnLayout, 'Text', 'Print Table', ...
                   'ButtonPushedFcn', @(btn,event) printTable(uitTable));
    printBtn.Layout.Row = 2;
    printBtn.Layout.Column = 2;

    statsBtn = uibutton(btnLayout, 'Text', 'Show Statistics', ...
                    'ButtonPushedFcn', @(btn,event) showStatistics(uitTable));
    statsBtn.Layout.Row = 2;
    statsBtn.Layout.Column = 3;
    %}

    % Update sort dropdown with column names
updateSortDropdown(sortByDropdown, uitTable);

    % Create callback to update row count when table changes
addlistener(uitTable, 'Data', 'PostSet', @(src, event) updateRowCount(rowCountLabel, uitTable));
drawnow;
fig.Visible = 'on';
end

function views = i_normalizeviews(T)
% Accept a single table (or empty) as well as a struct array of named views.
if isstruct(T) && isfield(T, 'Name') && isfield(T, 'Table') && ~isempty(T)
    views = T(:)';
else
    views = struct('Name', 'Table', 'Table', {T});
end
end

function view = i_currentview(fig)
view = fig.UserData.Views(fig.UserData.CurrentIndex);
end

function [data, headers] = i_viewcontent(view)
T = view.Table;
if isempty(T) && ~istable(T)
    % No table supplied at all; fall back to sample data. An empty table is
    % still shown as its own headers with zero rows.
    data = generateSampleData();
    headers = {'Name', 'Age', 'Height (cm)', 'Weight (kg)', 'BMI'};
else
    T = convertvars(T, @isstring, 'cellstr');
    headers = T.Properties.VariableNames;
    data = table2cell(T);
end
end

function name = i_defname(fig, defname)
% Qualify the default export filename with the view being displayed so that
% exporting several views does not keep proposing the same file.
name = defname;
if isempty(name) || numel(fig.UserData.Views) < 2, return; end
% Drop any trailing row count, e.g. "Up-regulated (10)" -> "Up_regulated"
suffix = regexprep(char(string(i_currentview(fig).Name)), '\s*\(\d+\)\s*$', '');
name = sprintf('%s_%s', char(string(name)), matlab.lang.makeValidName(suffix));
end

function switchView(dropdown, fig, tableObj, searchField, rowCountLabel, sortByDropdown)
% Show the table selected in the View dropdown.
idx = find(strcmp({fig.UserData.Views.Name}, dropdown.Value), 1);
if isempty(idx), return; end
fig.UserData.CurrentIndex = idx;
searchField.Value = '';
[data, headers] = i_viewcontent(i_currentview(fig));
tableObj.Data = data;
tableObj.ColumnName = headers;
updateRowCount(rowCountLabel, tableObj);
updateSortDropdown(sortByDropdown, tableObj);
end

function data = generateSampleData()
% Generate some sample data for the table
names = {'John Smith', 'Mary Johnson', 'Robert Williams', 'Susan Brown', ...
         'Michael Davis', 'Patricia Miller', 'James Wilson', 'Linda Moore', ...
         'David Taylor', 'Jennifer Anderson', 'William Thomas', 'Elizabeth Jackson'};
ages = randi([18, 65], length(names), 1);
heights = randi([150, 190], length(names), 1);
weights = randi([50, 100], length(names), 1);

% Calculate BMI (weight in kg / (height in m)^2)
bmi = weights ./ ((heights/100).^2);
bmi = round(bmi * 10) / 10;  % Round to 1 decimal place

% Combine all data into a table
data = table(names', ages, heights, weights, bmi, ...
            'VariableNames', {'Name', 'Age', 'Height', 'Weight', 'BMI'});
data = table2cell(data);
end

function exportToCSV(tableObj, defname)
% Function to export table data to a CSV file
[file, path] = uiputfile('*.csv', 'Save Table Data', defname);
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

function exportToExcel(tableObj, defname)
% Function to export table data to an Excel file
[file, path] = uiputfile({'*.xlsx;*.xls', 'Excel Files (*.xlsx, *.xls)'}, ...
'Save Table Data', defname);
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

function exportToMAT(tableObj, defname)
% Function to export table data to a MAT file
[file, path] = uiputfile('*.mat', 'Save Table Data', defname);
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

function exportToWorkspace(tableObj, fig)
T = i_tableFromUITable(tableObj);
labels = {'Save to variable named:'};
vars = {'T'};
values = {T};
gui.i_bringtofront(fig);
gui.myExport2wsdlg(labels, vars, values, ...
'Save Data to Workspace', [], fig);
end

function executeCustomAction(action, tableObj, fig)
try
    if nargin(action.Callback) >= 2
        action.Callback(tableObj, fig);
    else
        action.Callback();
    end
catch ME
    gui.myErrordlg(fig, ME.message, ME.identifier);
end
end

function T = i_tableFromUITable(tableObj)
columnNames = tableObj.ColumnName;
data = tableObj.Data;
T = cell2table(data, 'VariableNames', strrep(columnNames, ' ', '_'));
end

function refreshData(tableObj, rowCountLabel, sortByDropdown, fig)
% Function to reload the displayed view from its source table
[data, headers] = i_viewcontent(i_currentview(fig));
tableObj.Data = data;
tableObj.ColumnName = headers;

% Update row count
updateRowCount(rowCountLabel, tableObj);

% Reset sort dropdown
updateSortDropdown(sortByDropdown, tableObj);
end

%{
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
uiexportdlg(printFig);
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
%}

function searchTable(searchField, tableObj, fig)
% Function to search the displayed view
searchText = lower(searchField.Value);

% Always search the full view, not the currently filtered rows
data = i_viewcontent(i_currentview(fig));

% If search text is empty, show all rows
if isempty(searchText)
    tableObj.Data = data;
    return;
end

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
