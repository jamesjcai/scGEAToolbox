function callback_CellGeneProfile(src, ~)
%CALLBACK_CELLGENEPROFILE Colour every gene by how one clicked cell differs from the rest.
%
%   Behind Analyze > Gene Statistics > Gene Profile of a Clicked Cell. The
%   left panel is the embedding; clicking a cell makes it the focused cell.
%   The right panel draws one point per gene, mean expression against the
%   fraction of cells that detect it, coloured by how much higher or lower
%   the gene is in the focused cell than in the average cell. The cell is
%   averaged with its nearest neighbours in principal-component space
%   (PKG.E_CELLPROFILEROW), because most genes are not detected in any one
%   cell; the neighbours are marked on the embedding.
%
%   The subtitle names the genes most raised and most lowered. Toolbar
%   buttons switch the colour to the neighbourhood's own expression and
%   back, set how many cells are averaged, step back and forward through the
%   cells focused so far, find a cell by its ID, list every gene in a table,
%   and send the values to the workspace.
%
%   The idea is AnnZarro's gene plot coloured by the focused cell's row of
%   a cells x genes layer (Otto, Baasri and Setty, 2026).
%
%   See also PKG.E_CELLPROFILEMODEL, PKG.E_CELLPROFILEROW,
%   GUI.CALLBACK_CELLNEIGHBOURHOOD.

[FigureHandle, sce] = gui.gui_getfigsce(src);
if isempty(gui.i_checkexistingembed(sce))
    gui.myWarndlg(FigureHandle, ['The cells have no embedding yet. ' ...
        'Embed them first, then try again.']);
    return;
end
if sce.NumCells < 3 || sce.NumGenes < 2
    gui.myWarndlg(FigureHandle, 'At least three cells and two genes are needed.');
    return;
end

fw = gui.myWaitbar(FigureHandle);
try
    model = pkg.e_cellprofilemodel(sce.X, sce.g);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

numNamed = 5;
numLabelled = 5;
maxNeighbours = 500;
blue = [0.02, 0.19, 0.38];
white = [0.97, 0.97, 0.97];
red = [0.40, 0, 0.05];
diffMap = interp1([0; 0.5; 1], [blue; white; red], linspace(0, 1, 255)');
lightGrey = [0.88, 0.88, 0.88];
purple = [0.25, 0, 0.45];
levelMap = interp1([0; 1], [lightGrey; purple], linspace(0, 1, 256)');
neighbourColor = [0.13, 0.4, 0.67];

genes = model.Genes;
isShown = model.Fraction > 0;
xGene = model.Mean;
yGene = model.Fraction;
cellIds = string(sce.c_cell_id);
groupLabels = string(sce.c_cell_type_tx);
if numel(unique(groupLabels)) < 2
    groupLabels = string(sce.c_cluster_id);
end

s = sce.s;
is3d = size(s, 2) >= 3;
dotSize = 20;
if sce.NumCells > 5000
    dotSize = 8;
end

numNeighbours = min(15, sce.NumCells);
showDifference = true;
difference = [];
level = [];
history = 1;
historyPos = 1;

hx = gui.myFigure(FigureHandle);
hFig = hx.FigHandle;
hFig.Name = 'Gene Profile of a Clicked Cell';
hFig.Position(3) = max(hFig.Position(3), 1100);
showButton = hx.addCustomButton('off', @in_toggleshow, 'list.gif', ...
    'Colour by the neighbourhood''s own expression');
hx.addCustomButton('off', @in_setk, 'group_work_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', ...
    'Set how many cells are averaged...');
backIcon = gui.myFigure.backIcon();
hx.addCustomButton('on', @(~, ~) in_stephistory(-1), backIcon, ...
    'Previous focused cell');
hx.addCustomButton('off', @(~, ~) in_stephistory(1), fliplr(backIcon), ...
    'Next focused cell');
hx.addCustomButton('off', @in_findcell, 'barcode.jpg', ...
    'Focus a cell by its ID...');
hx.addCustomButton('on', @in_showtable, 'table.gif', ...
    'List every gene in a table');
hx.addCustomButton('off', @in_export, 'export.gif', ...
    'Send the values to the workspace...');

tl = tiledlayout(hFig, 1, 2, TileSpacing="compact", Padding="compact");
axCell = nexttile(tl);
if is3d
    hCells = scatter3(axCell, s(:, 1), s(:, 2), s(:, 3), dotSize, lightGrey, 'filled');
else
    hCells = scatter(axCell, s(:, 1), s(:, 2), dotSize, lightGrey, 'filled');
end
hCells.ButtonDownFcn = @in_clickcell;
hold(axCell, 'on');
hNeighbours = plot(axCell, nan, nan, '.', Color=neighbourColor, MarkerSize=12, ...
    PickableParts='none', HitTest='off');
hFocus = plot(axCell, nan, nan, 'o', MarkerSize=9, ...
    MarkerFaceColor=[0.85, 0.1, 0.1], MarkerEdgeColor='k', LineWidth=1.5, ...
    PickableParts='none', HitTest='off');
if is3d
    hNeighbours.ZData = nan;
    hFocus.ZData = nan;
end
hold(axCell, 'off');
box(axCell, 'on');
axis(axCell, 'tight');

axGene = nexttile(tl);
hGenes = scatter(axGene, xGene(isShown), yGene(isShown), 12, ...
    zeros(nnz(isShown), 1), 'filled');
hLabels = gobjects(0);
cbGene = colorbar(axGene);
xlabel(axGene, 'Mean expression over all cells (log1p, library-size normalised)');
ylabel(axGene, 'Fraction of cells detecting the gene');
box(axGene, 'on');

% Start on the cell nearest the middle of the embedding.
[~, history] = min(vecnorm(s - median(s), 2, 2));
in_draw();
hx.show(FigureHandle);

    function in_clickcell(~, evt)
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
        [difference, level, neighbours] = pkg.e_cellprofilerow(model, idx, ...
            numNeighbours);

        hNeighbours.XData = s(neighbours, 1);
        hNeighbours.YData = s(neighbours, 2);
        hFocus.XData = s(idx, 1);
        hFocus.YData = s(idx, 2);
        if is3d
            hNeighbours.ZData = s(neighbours, 3);
            hFocus.ZData = s(idx, 3);
        end
        title(axCell, sprintf('Focused cell %s (%s)', cellIds(idx), ...
            groupLabels(idx)), Interpreter='none');
        if numNeighbours > 1
            subtitle(axCell, sprintf('Averaged with its %d nearest cells (blue)', ...
                numNeighbours - 1));
        else
            subtitle(axCell, 'The cell alone, not averaged');
        end

        if showDifference
            value = difference;
            [~, order] = sort(abs(value(isShown)), 'ascend');
            colormap(axGene, diffMap);
            top = max(abs(value(isShown)));
            if ~(top > 0)
                top = 1;
            end
            clim(axGene, [-top, top]);
            cbGene.Label.String = 'Difference from the average cell (log1p)';
        else
            value = level;
            [~, order] = sort(value(isShown), 'ascend');
            colormap(axGene, levelMap);
            top = max(value(isShown));
            if ~(top > 0)
                top = 1;
            end
            clim(axGene, [0, top]);
            cbGene.Label.String = 'Expression in the neighbourhood (log1p)';
        end
        shown = find(isShown);
        shown = shown(order);
        hGenes.XData = xGene(shown);
        hGenes.YData = yGene(shown);
        hGenes.CData = value(shown);

        [~, upOrder] = sort(difference, 'descend');
        [~, downOrder] = sort(difference, 'ascend');
        up = upOrder(1:min(numNamed, numel(upOrder)));
        down = downOrder(1:min(numNamed, numel(downOrder)));
        delete(hLabels);
        labelled = [up(1:min(numLabelled, end)); down(1:min(numLabelled, end))];
        hLabels = text(axGene, xGene(labelled), yGene(labelled), ...
            "  " + genes(labelled), FontSize=8, Interpreter='none', ...
            PickableParts='none', HitTest='off');
        title(axGene, sprintf('Genes in %s against the average cell', ...
            cellIds(idx)), Interpreter='none');
        subtitle(axGene, {sprintf('Raised: %s', in_namelist(up)), ...
            sprintf('Lowered: %s', in_namelist(down))}, Interpreter='none');
    end

    function txt = in_namelist(list)
        txt = strjoin(compose('%s %+.2f', genes(list), difference(list)), ', ');
    end

    function in_toggleshow(~, ~)
        showDifference = ~showDifference;
        if showDifference
            showButton.TooltipString = 'Colour by the neighbourhood''s own expression';
        else
            showButton.TooltipString = 'Colour by the difference from the average cell';
        end
        in_draw();
    end

    function in_setk(~, ~)
        k = gui.i_inputnumk(numNeighbours, 1, min(maxNeighbours, sce.NumCells), ...
            'Number of cells to average, the focused cell included', hFig);
        if isempty(k)
            return;
        end
        numNeighbours = k;
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

    function T = in_table()
        T = table(genes, difference, level, model.Mean, model.Fraction, ...
            VariableNames=["Gene", "Difference", "NeighbourhoodExpression", ...
            "MeanExpression", "FractionDetected"]);
        T = sortrows(T, "Difference", "descend");
    end

    function in_showtable(~, ~)
        T = in_table();
        gui.i_openinstead(hFig, @() gui.TableViewerApp(T, hFig, ...
            "CellGeneProfile"));
    end

    function in_export(~, ~)
        idx = history(historyPos);
        labels = {'Gene table:', 'Focused cell index:'};
        vars = {'Tcellprofile', 'focusedcell'};
        gui.i_export2wsdlg(hFig, labels, vars, {in_table(), idx}, ...
            'Save Cell Gene Profile to Workspace');
    end
end
