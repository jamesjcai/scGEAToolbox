function [d, tf] = myExport2wsdlg(labels, vars, vals, titleText, defs, parentfig)
    % Validate input arguments
    tf = 1;
    narginchk(3, 6);
    if nargin < 6
        parentfig = [];
    end
    if nargin < 5 || isempty(defs)
        defs = true(1, length(labels));
    end    
    if nargin < 4
        titleText = 'Export Variables to Workspace';
    end
    assert(iscell(labels) && iscell(vars) && iscell(vals), ...
        'Labels, vars, and vals must be cell arrays.');
    assert(numel(labels) == numel(vars) && numel(vars) == numel(vals), ...
        'Labels, vars, and vals must have the same number of elements.');


     dialogWidth = 400;
     dialogHeight = 100 + 30 * numel(labels); % Adjust height based on the number of inputs

    % Compute center position relative to parentfig
    parentPos = parentfig.Position;
    dialogX = parentPos(1) + (parentPos(3) - dialogWidth) / 2;
    dialogY = parentPos(2) + (parentPos(4) - dialogHeight) / 2;
    % d = uifigure('Position', [dialogX, dialogY, dialogWidth, dialogHeight], ...
    %              'Name', titleText, 'WindowStyle', 'modal');

    % Create the dialog window
    d = uifigure('Name', titleText, 'WindowStyle', 'modal');
    d.Position = [dialogX, dialogY, dialogWidth, dialogHeight];

    % Create a grid layout manager
    gl = uigridlayout(d, [numel(labels) + 1, 3]);
    gl.RowSpacing = 5;
    gl.ColumnSpacing = 4;
    gl.Padding = [10, 10, 10, 10];
    gl.ColumnWidth = {'1x', '1x', 100}; % Corrected property name

    % Initialize components
    checkboxes = gobjects(numel(labels), 1);
    editfields = gobjects(numel(labels), 1);

    % Populate the grid with labels, checkboxes, and edit fields
    for i = 1:numel(labels)
        checkboxes(i) = uicheckbox(gl, 'Value', defs(i), "Text","");
        checkboxes(i).Layout.Row = i;
        checkboxes(i).Layout.Column = 4;

        lbl = uilabel(gl, 'Text', labels{i}, 'HorizontalAlignment', 'right');
        lbl.Layout.Row = i;
        lbl.Layout.Column = [1, 2];

        a = pkg.ai_suggetVarName(vars{i}, 'base');
        

        editfields(i) = uieditfield(gl, 'text', 'Value', a);
        editfields(i).Layout.Row = i;
        editfields(i).Layout.Column = 3;
    end

    % Add OK and Cancel buttons
    btnOK = uibutton(gl, 'Text', 'OK', 'ButtonPushedFcn', @(~, ~) onOK());
    btnOK.Layout.Row = numel(labels) + 1;
    btnOK.Layout.Column = 2;

    btnCancel = uibutton(gl, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) delete(d));
    btnCancel.Layout.Row = numel(labels) + 1;
    btnCancel.Layout.Column = 3;

    % Wait for the dialog to close before returning
    uiwait(d);

    function onOK()
        % Export selected variables to the base workspace
        for idx = 1:numel(labels)
            if checkboxes(idx).Value
                varName = editfields(idx).Value;
                if isvarname(varName)
                    assignin('base', varName, vals{idx});
                else
                    uialert(d, ['Invalid variable name: ', varName], 'Error');
                    return;
                end
            end
        end
        delete(d);
    end
end
