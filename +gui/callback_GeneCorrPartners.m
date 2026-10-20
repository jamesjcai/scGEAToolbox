function callback_GeneCorrPartners(src, ~)
%CALLBACK_GENECORRPARTNERS Colour every gene by its correlation with one clicked gene.
%
%   Behind Analyze > Gene Statistics > Gene Correlation Partners. The left
%   panel draws one point per gene, mean expression against the fraction of
%   cells that detect it, coloured by the gene's correlation across cells
%   with the focused gene, on a fixed -1 to 1 scale so that one colour means
%   one correlation from gene to gene. Clicking a gene makes it the focused
%   gene. The right panel colours the embedding by the focused gene's
%   expression, so its partners can be read against where it is expressed.
%
%   The subtitle names the strongest partners. Toolbar buttons find a gene
%   by name, step back and forward through the genes focused so far, switch
%   between Pearson (log1p normalised expression) and Spearman correlation,
%   fit the colour range to the data, list every gene's correlation in a
%   table, and send the correlations to the workspace.
%
%   The idea is AnnZarro's focus model (Otto, Baasri and Setty, 2026), where
%   a focused gene picks a row of a genes x genes matrix. Here the row is
%   computed when the gene is clicked, so nothing genes x genes is stored.
%
%   See also PKG.E_GENECORRMODEL, PKG.E_GENECORRROW,
%   GUI.CALLBACK_CELLNEIGHBOURHOOD.

[FigureHandle, sce] = gui.gui_getfigsce(src);
if sce.NumGenes < 2 || sce.NumCells < 3
    gui.myWarndlg(FigureHandle, 'At least two genes and three cells are needed.');
    return;
end

fw = gui.myWaitbar(FigureHandle);
try
    model = pkg.e_genecorrmodel(sce.X, sce.g);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

numPartnersNamed = 5;
numPartnersLabelled = 5;
minFractionToStart = 0.05;
blue = [0.02, 0.19, 0.38];
white = [0.97, 0.97, 0.97];
red = [0.40, 0, 0.05];
corrMap = interp1([0; 0.5; 1], [blue; white; red], linspace(0, 1, 255)');
lightGrey = [0.92, 0.92, 0.92];
purple = [0.25, 0, 0.45];
exprMap = interp1([0; 1], [lightGrey; purple], linspace(0, 1, 256)');

genes = model.Genes;
isShown = model.Fraction > 0;
xGene = model.Mean;
yGene = model.Fraction;

hasEmbedding = ~isempty(gui.i_checkexistingembed(sce));
s = sce.s;
is3d = size(s, 2) >= 3;
dotSize = 20;
if sce.NumCells > 5000
    dotSize = 8;
end

method = "pearson";
isRangeFitted = false;
rho = [];
history = 1;
historyPos = 1;

hx = gui.myFigure(FigureHandle);
hFig = hx.FigHandle;
hFig.Name = 'Gene Correlation Partners';
hFig.Position(3) = max(hFig.Position(3), 1100);
hx.addCustomButton('off', @in_findgene, 'label.gif', ...
    'Focus a gene by its name...');
backIcon = gui.myFigure.backIcon();
hx.addCustomButton('off', @(~, ~) in_stephistory(-1), backIcon, ...
    'Previous focused gene');
hx.addCustomButton('off', @(~, ~) in_stephistory(1), fliplr(backIcon), ...
    'Next focused gene');
methodButton = hx.addCustomButton('on', @in_togglemethod, 'list.gif', ...
    'Switch to Spearman correlation');
rangeButton = hx.addCustomButton('off', @in_togglerange, 'setrange.gif', ...
    'Fit the colour range to the data');
hx.addCustomButton('on', @in_showtable, 'table.gif', ...
    'List every gene''s correlation in a table');
hx.addCustomButton('off', @in_export, 'export.gif', ...
    'Send the correlations to the workspace...');

tl = tiledlayout(hFig, 1, 1 + hasEmbedding, TileSpacing="compact", ...
    Padding="compact");
axGene = nexttile(tl);
hGenes = scatter(axGene, xGene(isShown), yGene(isShown), 12, ...
    zeros(nnz(isShown), 1), 'filled');
hGenes.ButtonDownFcn = @in_clickgene;
hold(axGene, 'on');
hFocus = plot(axGene, nan, nan, 'o', MarkerSize=10, MarkerEdgeColor='k', ...
    LineWidth=1.5, PickableParts='none', HitTest='off');
hold(axGene, 'off');
hLabels = gobjects(0);
colormap(axGene, corrMap);
cbGene = colorbar(axGene);
xlabel(axGene, 'Mean expression (log1p, library-size normalised)');
ylabel(axGene, 'Fraction of cells detecting the gene');
box(axGene, 'on');

if hasEmbedding
    axCell = nexttile(tl);
    if is3d
        hCells = scatter3(axCell, s(:, 1), s(:, 2), s(:, 3), dotSize, ...
            zeros(sce.NumCells, 1), 'filled');
    else
        hCells = scatter(axCell, s(:, 1), s(:, 2), dotSize, ...
            zeros(sce.NumCells, 1), 'filled');
    end
    colormap(axCell, exprMap);
    cbCell = colorbar(axCell);
    cbCell.Label.String = 'log1p normalised expression';
    box(axCell, 'on');
    axis(axCell, 'tight');
end

% Start on the most variable gene among those most cells detect.
spread = model.Pearson.Std;
spread(model.Fraction < minFractionToStart) = -Inf;
[~, history] = max(spread);
in_draw();
hx.show(FigureHandle);

    function in_clickgene(~, evt)
        if ~strcmp(hFig.SelectionType, 'normal')
            return;
        end
        pt = evt.IntersectionPoint;
        % Nearest on screen, so each axis is scaled by its range.
        dx = (xGene - pt(1))/diff(axGene.XLim);
        dy = (yGene - pt(2))/diff(axGene.YLim);
        d = dx.^2 + dy.^2;
        d(~isShown) = Inf;
        [~, idx] = min(d);
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
        rho = pkg.e_genecorrrow(model, idx, method);

        % Strong correlations of either sign are drawn last, on top.
        shown = find(isShown & ~isnan(rho));
        strength = abs(rho(shown));
        strength(shown == idx) = Inf;
        [~, order] = sort(strength, 'ascend');
        shown = shown(order);
        hGenes.XData = xGene(shown);
        hGenes.YData = yGene(shown);
        hGenes.CData = rho(shown);
        hFocus.XData = xGene(idx);
        hFocus.YData = yGene(idx);

        if isRangeFitted
            others = rho;
            others(idx) = NaN;
            top = max(abs(others), [], 'omitnan');
            if ~(top > 0)
                top = 1;
            end
            clim(axGene, [-top, top]);
        else
            clim(axGene, [-1, 1]);
        end
        cbGene.Label.String = sprintf('%s correlation with %s', ...
            in_methodname(), genes(idx));

        [partners, partnerRho] = in_partners(idx);
        delete(hLabels);
        numLabels = min(numPartnersLabelled, numel(partners));
        hLabels = text(axGene, xGene(partners(1:numLabels)), ...
            yGene(partners(1:numLabels)), "  " + genes(partners(1:numLabels)), ...
            FontSize=8, Interpreter='none', PickableParts='none', HitTest='off');
        title(axGene, sprintf('Genes correlated with %s', genes(idx)), ...
            Interpreter='none');
        numNamed = min(numPartnersNamed, numel(partners));
        named = compose('%s %.2f', genes(partners(1:numNamed)), ...
            partnerRho(1:numNamed));
        subtitle(axGene, sprintf('Strongest partners: %s', ...
            strjoin(named, ', ')), Interpreter='none');

        if hasEmbedding
            expr = full(model.Pearson.Matrix(idx, :)).';
            [~, cellOrder] = sort(expr, 'ascend');
            hCells.XData = s(cellOrder, 1);
            hCells.YData = s(cellOrder, 2);
            if is3d
                hCells.ZData = s(cellOrder, 3);
            end
            hCells.CData = expr(cellOrder);
            top = max(expr);
            if ~(top > 0)
                top = 1;
            end
            clim(axCell, [0, top]);
            title(axCell, sprintf('%s expression', genes(idx)), Interpreter='none');
            subtitle(axCell, sprintf('Detected in %.1f%% of cells', ...
                100*model.Fraction(idx)));
        end
    end

    function [partners, partnerRho] = in_partners(idx)
        % Other genes by correlation, strongest positive first.
        candidates = rho;
        candidates(idx) = NaN;
        [partnerRho, partners] = sort(candidates, 'descend', ...
            MissingPlacement='last');
        keep = ~isnan(partnerRho);
        partners = partners(keep);
        partnerRho = partnerRho(keep);
    end

    function name = in_methodname()
        if method == "spearman"
            name = 'Spearman';
        else
            name = 'Pearson';
        end
    end

    function in_findgene(~, ~)
        answer = gui.myInputdlg({'Gene name:'}, 'Find Gene', ...
            {char(genes(history(historyPos)))}, hFig);
        if isempty(answer)
            return;
        end
        query = strtrim(string(answer{1}));
        idx = find(strcmpi(genes, query), 1);
        if isempty(idx)
            idx = find(startsWith(genes, query, IgnoreCase=true), 1);
        end
        if isempty(idx)
            gui.myWarndlg(hFig, sprintf('No gene name starts with "%s".', query));
            return;
        end
        in_setfocus(idx);
    end

    function in_togglemethod(~, ~)
        if method == "pearson"
            if isempty(model.Spearman)
                fwRank = gui.myWaitbar(hFig);
                try
                    model = pkg.e_genecorrmodel(sce.X, sce.g, Spearman=true);
                catch ME
                    gui.myWaitbar(hFig, fwRank, true);
                    gui.myErrordlg(hFig, ME.message, ME.identifier);
                    return;
                end
                gui.myWaitbar(hFig, fwRank);
            end
            method = "spearman";
            methodButton.TooltipString = 'Switch to Pearson correlation';
        else
            method = "pearson";
            methodButton.TooltipString = 'Switch to Spearman correlation';
        end
        in_draw();
    end

    function in_togglerange(~, ~)
        isRangeFitted = ~isRangeFitted;
        if isRangeFitted
            rangeButton.TooltipString = 'Fix the colour range at -1 to 1';
        else
            rangeButton.TooltipString = 'Fit the colour range to the data';
        end
        in_draw();
    end

    function T = in_table()
        idx = history(historyPos);
        T = table(genes, rho, model.Mean, model.Fraction, ...
            VariableNames=["Gene", "Correlation", "MeanExpression", ...
            "FractionDetected"]);
        T(idx, :) = [];
        T = sortrows(T, "Correlation", "descend", MissingPlacement="last");
    end

    function in_showtable(~, ~)
        T = in_table();
        gui.i_openinstead(hFig, @() gui.TableViewerApp(T, hFig, ...
            "GeneCorrelationPartners"));
    end

    function in_export(~, ~)
        labels = {'Correlation table:', 'Correlation with each gene, in gene order:'};
        vars = {'Tgenecorr', 'genecorr'};
        gui.i_export2wsdlg(hFig, labels, vars, {in_table(), rho}, ...
            'Save Gene Correlations to Workspace');
    end
end
