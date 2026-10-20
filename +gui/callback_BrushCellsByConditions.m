function callback_BrushCellsByConditions(src, ~)
%CALLBACK_BRUSHCELLSBYCONDITIONS Pick cells by conditions on their attributes.
%
%   Behind Edit > Brush Cells by Conditions. A window lists conditions
%   on per-cell variables - cell type, cluster, batch, QC metrics and every
%   cell attribute (PKG.E_CELLCONDITIONTABLE) - such as "CellType Equals
%   Beta cells" and "MitoPercent Less than 5", met all together (AND) or
%   any one of them (OR). As the conditions change, the embedding beside
%   them keeps the matching cells in their group colours and greys out the
%   rest, or hides them, and a count says how many match.
%
%   Brush in Main Window makes the matching cells the brushed selection
%   of the main plot, so every command that acts on brushed cells - find
%   marker genes, add them to a new group, delete them - acts on them.
%   Export Matching Cells saves their rows of the table. Nothing is removed
%   from the data here.
%
%   The idea is AnnZarro's table filter (Otto, Baasri and Setty, 2026),
%   where a filtered cell table greys out the cells it leaves out in a
%   linked plot.
%
%   See also PKG.E_MATCHCONDITIONS, GUI.I_SETBRUSHDATA,
%   GUI.CALLBACK_SELECTCELLSBYCLASS.

[FigureHandle, sce] = gui.gui_getfigsce(src);
if sce.NumCells == 0
    gui.myWarndlg(FigureHandle, 'No data loaded.');
    return;
end

T = pkg.e_cellconditiontable(sce);
varNames = string(T.Properties.VariableNames);
isNumericVar = varfun(@isnumeric, T, OutputFormat="uniform");
maxListedValues = 200;
rowHeight = 26;
lightGrey = [0.85, 0.85, 0.85];

numCells = sce.NumCells;
hasEmbedding = ~isempty(gui.i_checkexistingembed(sce));
s = sce.s;
is3d = size(s, 2) >= 3;
dotSize = 20;
if numCells > 5000
    dotSize = 8;
end
[groupId, ~] = findgroups(string(sce.c(:)));
groupColors = pkg.i_mycolorlines(max(groupId));

conditions = struct('Variable', {}, 'Operator', {}, 'Value', {});
combine = "all";
hideOthers = false;
mask = true(numCells, 1);
% False while a condition cannot be evaluated; MASK is then stale.
isValid = true;

fig = uifigure(Name='Brush Cells by Conditions', ...
    Position=[0, 0, 1050, 560], Visible='off');
gui.i_movegui2parent(fig, FigureHandle);
mainGrid = uigridlayout(fig, [1, 2], ColumnWidth={470, '1x'});
leftGrid = uigridlayout(mainGrid, [6, 1], ...
    RowHeight={rowHeight, '1x', 30, 40, rowHeight, 30});
uidropdown(leftGrid, ...
    Items={'Show cells that meet all conditions (AND)', ...
    'Show cells that meet any condition (OR)'}, ...
    ItemsData={'all', 'any'}, ValueChangedFcn=@in_setcombine);
condPanel = uipanel(leftGrid, Title='Conditions');
condGrid = uigridlayout(condPanel, [1, 4], ...
    ColumnWidth={'1x', 110, '1x', 26}, Scrollable='on');
uibutton(leftGrid, Text='Add Condition', ButtonPushedFcn=@in_add);
countLabel = uilabel(leftGrid, WordWrap='on', FontWeight='bold');
uicheckbox(leftGrid, Text='Hide cells that do not match', ...
    ValueChangedFcn=@in_sethide);
buttonGrid = uigridlayout(leftGrid, [1, 2], Padding=[0, 0, 0, 0]);
uibutton(buttonGrid, Text='Brush in Main Window', ...
    Tooltip='Make the matching cells the brushed selection of the main plot', ...
    ButtonPushedFcn=@in_brush);
uibutton(buttonGrid, Text='Export Matching Cells...', ...
    ButtonPushedFcn=@in_export);

if hasEmbedding
    ax = uiaxes(mainGrid);
    if is3d
        hOthers = scatter3(ax, nan, nan, nan, dotSize, lightGrey, 'filled');
        hold(ax, 'on');
        hMatch = scatter3(ax, nan, nan, nan, dotSize, [0, 0, 0], 'filled');
        view(ax, 3);
    else
        hOthers = scatter(ax, nan, nan, dotSize, lightGrey, 'filled');
        hold(ax, 'on');
        hMatch = scatter(ax, nan, nan, dotSize, [0, 0, 0], 'filled');
    end
    hold(ax, 'off');
    box(ax, 'on');
    ax.XLim = in_limits(s(:, 1));
    ax.YLim = in_limits(s(:, 2));
    if is3d
        ax.ZLim = in_limits(s(:, 3));
    end
else
    uilabel(mainGrid, Text='The cells have no embedding to show.', ...
        HorizontalAlignment='center');
end

in_add();
fig.Visible = 'on';

    function lim = in_limits(x)
        lim = [min(x), max(x)];
        if ~(lim(2) > lim(1))
            lim = lim + [-1, 1];
        end
    end

    function in_add(~, ~)
        if isempty(conditions)
            name = varNames(1);
            if ismember("CellType", varNames)
                name = "CellType";
            end
        else
            name = string(conditions(end).Variable);
        end
        [op, value] = in_defaults(name);
        conditions(end + 1) = struct('Variable', name, 'Operator', op, ...
            'Value', value);
        in_rebuild();
        in_update();
    end

    function [op, value] = in_defaults(name)
        x = T.(name);
        if isNumericVar(varNames == name)
            op = "Greater than";
            value = string(sprintf('%.4g', median(x, 'omitnan')));
        else
            op = "Equals";
            values = in_listedvalues(name);
            value = values(1);
        end
    end

    function values = in_listedvalues(name)
        % The text values of a column, most frequent first.
        x = T.(name);
        x = x(~ismissing(x));
        [u, ~, ic] = unique(x);
        if isempty(u)
            values = "";
            return;
        end
        counts = accumarray(ic, 1);
        [~, order] = sort(counts, 'descend');
        values = u(order(1:min(maxListedValues, numel(order))));
    end

    function in_rebuild()
        delete(condGrid.Children);
        numConditions = numel(conditions);
        condGrid.RowHeight = repmat({rowHeight}, 1, max(numConditions, 1));
        if numConditions == 0
            note = uilabel(condGrid, Text=['No conditions, so every cell ' ...
                'is shown. Add one below.']);
            note.Layout.Row = 1;
            note.Layout.Column = [1, 4];
            return;
        end
        for j = 1:numConditions
            in_buildrow(j);
        end
    end

    function in_buildrow(j)
        name = string(conditions(j).Variable);
        isNumeric = isNumericVar(varNames == name);
        if isNumeric
            ops = pkg.e_conditionoperators("numeric");
        else
            ops = pkg.e_conditionoperators("text");
        end
        varDropdown = uidropdown(condGrid, Items=cellstr(varNames), ...
            Value=char(name), ValueChangedFcn=@(o, ~) in_setvariable(j, o.Value));
        opDropdown = uidropdown(condGrid, Items=cellstr(ops), ...
            Value=char(conditions(j).Operator), ...
            ValueChangedFcn=@(o, ~) in_setfield(j, 'Operator', o.Value));
        if isNumeric
            x = T.(name);
            valueControl = uieditfield(condGrid, 'text', ...
                Value=char(conditions(j).Value), ...
                Placeholder=sprintf('%.4g to %.4g', min(x), max(x)), ...
                Tooltip=sprintf('A number; the cells range from %.4g to %.4g', ...
                min(x), max(x)), ...
                ValueChangedFcn=@(o, ~) in_setfield(j, 'Value', o.Value));
        else
            valueControl = uidropdown(condGrid, Editable='on', ...
                Items=cellstr(in_listedvalues(name)), ...
                Value=char(conditions(j).Value), ...
                Tooltip='Pick a value, or type any text', ...
                ValueChangedFcn=@(o, ~) in_setfield(j, 'Value', o.Value));
        end
        removeButton = uibutton(condGrid, Text=char(8722), ...
            Tooltip='Remove this condition', ...
            ButtonPushedFcn=@(~, ~) in_remove(j));
        controls = {varDropdown, opDropdown, valueControl, removeButton};
        for col = 1:numel(controls)
            controls{col}.Layout.Row = j;
            controls{col}.Layout.Column = col;
        end
    end

    function in_setvariable(j, name)
        name = string(name);
        [op, value] = in_defaults(name);
        conditions(j) = struct('Variable', name, 'Operator', op, 'Value', value);
        in_rebuild();
        in_update();
    end

    function in_setfield(j, field, value)
        conditions(j).(field) = string(value);
        in_update();
    end

    function in_remove(j)
        conditions(j) = [];
        in_rebuild();
        in_update();
    end

    function in_setcombine(o, ~)
        combine = string(o.Value);
        in_update();
    end

    function in_sethide(o, ~)
        hideOthers = o.Value;
        in_draw();
    end

    function in_update()
        try
            mask = pkg.e_matchconditions(T, conditions, combine);
        catch ME
            isValid = false;
            countLabel.Text = ME.message;
            countLabel.FontColor = [0.8, 0, 0];
            in_draw();
            return;
        end
        isValid = true;
        countLabel.FontColor = [0, 0, 0];
        countLabel.Text = sprintf('%d of %d cells match (%.1f%%).', ...
            nnz(mask), numCells, 100*nnz(mask)/numCells);
        in_draw();
    end

    function in_draw()
        if ~hasEmbedding
            return;
        end
        shownMatch = mask & isValid;
        others = ~shownMatch & ~(hideOthers && isValid);
        hOthers.XData = s(others, 1);
        hOthers.YData = s(others, 2);
        hMatch.XData = s(shownMatch, 1);
        hMatch.YData = s(shownMatch, 2);
        if is3d
            hOthers.ZData = s(others, 3);
            hMatch.ZData = s(shownMatch, 3);
        end
        hMatch.CData = groupColors(groupId(shownMatch), :);
        if isValid
            title(ax, sprintf('%d of %d cells match', nnz(mask), numCells));
        else
            title(ax, 'Fix the condition in red');
        end
    end

    function in_brush(~, ~)
        if ~in_checkvalid('Brush Cells')
            return;
        end
        h = [];
        if isa(src, 'matlab.apps.AppBase') && isprop(src, 'h')
            h = src.h;
        elseif pkg.i_isvalid(FigureHandle)
            h = findobj(FigureHandle, 'Type', 'Scatter');
        end
        if numel(h) > 1
            h = h(1);
        end
        if ~gui.i_setbrushdata(h, mask)
            uialert(fig, ['The main window does not show these cells: its ' ...
                'data may have changed since this window opened. Close ' ...
                'this window and open it again.'], 'Brush Cells');
            return;
        end
        figure(FigureHandle);
    end

    function in_export(~, ~)
        if ~in_checkvalid('Export Matching Cells')
            return;
        end
        if ~any(mask)
            uialert(fig, 'No cells match, so there is nothing to export.', ...
                'Export Matching Cells', Icon='info');
            return;
        end
        gui.i_exporttable(T(mask, :), true, 'Tcells', 'MatchingCells', ...
            [], [], fig);
    end

    function tf = in_checkvalid(dlgTitle)
        tf = isValid;
        if ~tf
            uialert(fig, ['A condition cannot be checked; its message is ' ...
                'shown in red. Fix or remove it first.'], dlgTitle);
        end
    end
end
