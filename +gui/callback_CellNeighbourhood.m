function callback_CellNeighbourhood(src, ~)
%CALLBACK_CELLNEIGHBOURHOOD Colour the embedding by nearness to one clicked cell.
%
%   Behind View > Show Cell Neighbourhood. A window draws the embedding;
%   clicking a cell makes it the focused cell, and every cell is recoloured
%   by how near it is to that one. Dark means near for every measure:
%
%     Diffusion walk       where a short random walk on the kNN graph from
%                          the focused cell ends up (PKG.E_CELLNEIGHBOURROW)
%     Expression distance  distance in principal-component space
%     Embedding distance   distance in the embedding drawn
%
%   The subtitle says where the walk's mass lands, or which groups the 30
%   nearest cells belong to. Toolbar buttons choose the measure and the walk
%   length, step back and forward through the cells focused so far, find a
%   cell by its ID, hold the colour range fixed so that one colour means one
%   value from cell to cell, and send the current values to the workspace.
%
%   The idea is AnnZarro's focus model (Otto, Baasri and Setty, 2026), where
%   a focused cell picks a row of a cells x cells matrix. Here the row is
%   computed when the cell is clicked, so nothing cells x cells is stored.
%
%   See also PKG.E_CELLNEIGHBOURMODEL, PKG.E_CELLNEIGHBOURROW.

[FigureHandle, sce] = gui.gui_getfigsce(src);
if isempty(gui.i_checkexistingembed(sce))
    gui.myWarndlg(FigureHandle, ['The cells have no embedding yet. ' ...
        'Embed them first, then try again.']);
    return;
end
if sce.NumCells < 3
    gui.myWarndlg(FigureHandle, 'At least three cells are needed.');
    return;
end

fw = gui.myWaitbar(FigureHandle);
try
    model = pkg.e_cellneighbourmodel(sce.X, sce.g, sce.s);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

metricNames = ["walk", "expression", "embedding"];
metricLabels = ["Diffusion walk", "Expression distance (PCA)", ...
    "Embedding distance"];
numNearest = 30;
lightGrey = [0.92, 0.92, 0.92];
darkBlue = [0.03, 0.19, 0.42];
nearMap = interp1([0; 1], [lightGrey; darkBlue], linspace(0, 1, 256)');

cellIds = string(sce.c_cell_id);
groupLabels = string(sce.c_cell_type_tx);
if numel(unique(groupLabels)) < 2
    groupLabels = string(sce.c_cluster_id);
end
[groupNames, ~, groupId] = unique(groupLabels);
hasGroups = numel(groupNames) > 1;

s = sce.s;
is3d = size(s, 2) >= 3;
dotSize = 20;
if sce.NumCells > 5000
    dotSize = 8;
end

metric = metricNames(1);
steps = 5;
isRangeLocked = false;
lockedRange = [];
values = [];
history = 1;
historyPos = 1;

hx = gui.myFigure(FigureHandle);
hFig = hx.FigHandle;
hFig.Name = 'Cell Neighbourhood';
hAx = hx.AxHandle;
hx.addCustomButton('off', @in_pickmetric, 'list.gif', ...
    'Choose how nearness is measured');
hx.addCustomButton('off', @in_setsteps, ...
    'turn_sharp_right_17dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', ...
    'Set the number of walk steps');
backIcon = gui.myFigure.backIcon();
hx.addCustomButton('on', @(~, ~) in_stephistory(-1), backIcon, ...
    'Previous focused cell');
hx.addCustomButton('off', @(~, ~) in_stephistory(1), fliplr(backIcon), ...
    'Next focused cell');
hx.addCustomButton('off', @in_findcell, 'barcode.jpg', ...
    'Focus a cell by its ID...');
lockButton = hx.addCustomButton('on', @in_togglelock, 'setrange.gif', ...
    'Lock colour range');
hx.addCustomButton('off', @in_export, 'export.gif', ...
    'Send the values to the workspace...');

if is3d
    hScatter = scatter3(hAx, s(:, 1), s(:, 2), s(:, 3), dotSize, ...
        zeros(sce.NumCells, 1), 'filled');
else
    hScatter = scatter(hAx, s(:, 1), s(:, 2), dotSize, ...
        zeros(sce.NumCells, 1), 'filled');
end
hScatter.ButtonDownFcn = @in_click;
hold(hAx, 'on');
hFocus = plot(hAx, nan, nan, 'o', MarkerSize=9, ...
    MarkerFaceColor=[0.85, 0.1, 0.1], MarkerEdgeColor='k', LineWidth=1.5, ...
    PickableParts='none', HitTest='off');
if is3d
    hFocus.ZData = nan;
end
hold(hAx, 'off');
hCbar = colorbar(hAx);
box(hAx, 'on');
axis(hAx, 'tight');
grid(hAx, 'off');

% Start on the cell nearest the middle of the embedding.
[~, history] = min(vecnorm(s - median(s), 2, 2));
in_draw();
hx.show(FigureHandle);

    function in_click(~, evt)
        if ~strcmp(hFig.SelectionType, 'normal')
            return;
        end
        pt = evt.IntersectionPoint;
        nd = 2 + is3d;
        [~, idx] = min(sum((s(:, 1:nd) - pt(1:nd)).^2, 2));
        in_setfocus(idx);
    end

    function in_setfocus(idx)
        if idx == history(historyPos)
            return;
        end
        history = [history(1:historyPos), idx];
        historyPos = numel(history);
        in_draw();
    end

    function in_stephistory(delta)
        newPos = historyPos + delta;
        if newPos < 1 || newPos > numel(history)
            return;
        end
        historyPos = newPos;
        in_draw();
    end

    function in_draw()
        idx = history(historyPos);
        values = pkg.e_cellneighbourrow(model, idx, metric, steps);
        isWalk = metric == "walk";

        % Near cells are drawn last, so the dense far majority does not
        % hide them.
        if isWalk
            [~, order] = sort(values, 'ascend');
        else
            [~, order] = sort(values, 'descend');
        end
        hScatter.XData = s(order, 1);
        hScatter.YData = s(order, 2);
        if is3d
            hScatter.ZData = s(order, 3);
        end
        hScatter.CData = values(order);
        hFocus.XData = s(idx, 1);
        hFocus.YData = s(idx, 2);
        if is3d
            hFocus.ZData = s(idx, 3);
        end

        if isWalk
            colormap(hAx, nearMap);
        else
            colormap(hAx, flipud(nearMap));
        end
        if isRangeLocked
            clim(hAx, lockedRange);
        else
            in_autorange(isWalk);
        end
        label = metricLabels(metricNames == metric);
        if isWalk
            label = sprintf('%s, %d steps', label, steps);
        end
        hCbar.Label.String = label;

        title(hAx, sprintf('Focused cell %s (%s)', cellIds(idx), ...
            groupLabels(idx)), Interpreter='none');
        subtitle(hAx, in_summary(idx, isWalk), Interpreter='none');
    end

    function in_autorange(isWalk)
        others = values;
        others(history(historyPos)) = [];
        top = max(others);
        bottom = min(others);
        if isWalk
            bottom = 0;
        end
        if ~(top > bottom)
            top = bottom + 1;
        end
        clim(hAx, [bottom, top]);
    end

    function txt = in_summary(idx, isWalk)
        if ~hasGroups
            txt = '';
            return;
        end
        if isWalk
            share = accumarray(groupId, values, [numel(groupNames), 1]);
            lead = 'Walk mass';
        else
            [~, nearOrder] = sort(values, 'ascend');
            nearOrder(nearOrder == idx) = [];
            near = nearOrder(1:min(numNearest, numel(nearOrder)));
            share = accumarray(groupId(near), 1, [numel(groupNames), 1]) / numel(near);
            lead = sprintf('Nearest %d cells', numel(near));
        end
        [share, top] = sort(share, 'descend');
        % Shares that would print as 0.0% say nothing.
        numShown = max(1, min(3, nnz(share >= 0.0005)));
        parts = compose('%.1f%% %s', 100*share(1:numShown), groupNames(top(1:numShown)));
        txt = sprintf('%s: %s', lead, strjoin(parts, ', '));
    end

    function in_pickmetric(~, ~)
        current = find(metricNames == metric);
        [indx, tf] = gui.myListdlg(hFig, cellstr(metricLabels), 'Measure', ...
            current, false, true, [300, 200], ...
            'How should nearness to the focused cell be measured?');
        if ~tf || isempty(indx) || indx == current
            return;
        end
        metric = metricNames(indx);
        % A range set on one measure means nothing on another.
        in_setlock(false);
        in_draw();
    end

    function in_setsteps(~, ~)
        k = gui.i_inputnumk(steps, 1, 500, 'Number of walk steps', hFig);
        if isempty(k)
            return;
        end
        steps = k;
        if metric ~= "walk"
            metric = "walk";
            in_setlock(false);
        end
        in_draw();
    end

    function in_findcell(~, ~)
        answer = gui.myInputdlg({'Cell ID (or the start of one):'}, ...
            'Find Cell', {char(cellIds(history(historyPos)))}, hFig);
        if isempty(answer)
            return;
        end
        query = strtrim(string(answer{1}));
        idx = find(cellIds == query, 1);
        if isempty(idx)
            idx = find(startsWith(cellIds, query, IgnoreCase=true), 1);
        end
        if isempty(idx)
            gui.myWarndlg(hFig, sprintf('No cell ID starts with "%s".', query));
            return;
        end
        in_setfocus(idx);
    end

    function in_togglelock(~, ~)
        in_setlock(~isRangeLocked);
        if ~isRangeLocked
            in_draw();
        end
    end

    function in_setlock(tf)
        isRangeLocked = tf;
        if tf
            lockedRange = clim(hAx);
            lockButton.TooltipString = sprintf(['Unlock colour range ' ...
                '(locked at %.4g to %.4g)'], lockedRange);
        else
            lockedRange = [];
            lockButton.TooltipString = 'Lock colour range';
        end
    end

    function in_export(~, ~)
        idx = history(historyPos);
        labels = {'Values, one per cell:', 'Focused cell index:'};
        vars = {'nearness', 'focusedcell'};
        gui.i_export2wsdlg(hFig, labels, vars, {values, idx}, ...
            'Save Cell Neighbourhood to Workspace');
    end
end
