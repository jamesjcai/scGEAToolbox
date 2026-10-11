function callback_CellStateDriverGenes(src, ~)
%CALLBACK_CELLSTATEDRIVERGENES Rank genes that drive cell-state structure.
%
%   Behind Analyze > Gene Statistics > Genes Driving Cell-State Structure.
%   Runs SC_POWERSEEK (Power-Seek, Ghosh et al., Cell Systems 2026) on the
%   brushed cells, one cell type, or all cells, and shows the cell-cell
%   covariance spectrum against its shuffled null, the choice of the
%   power-law regime, and the ranked genes.
%
%   The genes it ranks co-vary with the dominant structure among the cells.
%   The paper calls them memory genes, inherited across divisions; the same
%   signal also comes from lineage, a persistent transcriptional state or a
%   shared microenvironment, which one snapshot cannot tell apart.
%
%   See also SC_POWERSEEK

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('Power-Seek [doi:10.1016/j.cels.2026.101707]', ...
        FigureHandle)
    return;
end

[cellIdx, cellsName] = in_selectcells(src, sce, FigureHandle);
if isempty(cellIdx), return; end

defaultNumGenes = 3000;
answer = gui.myInputdlg({'Number of most variable genes to rank:', ...
    'Eigenvalues in the power-law regime (blank = choose automatically):'}, ...
    'Power-Seek', {num2str(defaultNumGenes), ''}, FigureHandle);
if isempty(answer), return; end
numGenes = str2double(answer{1});
if ~isfinite(numGenes) || numGenes < 10 || numGenes ~= round(numGenes)
    gui.myErrordlg(FigureHandle, 'The number of genes must be a whole number of at least 10.');
    return;
end
numEig = NaN;
if strlength(strtrim(string(answer{2}))) > 0
    numEig = str2double(answer{2});
    if ~isfinite(numEig) || numEig < 2 || numEig ~= round(numEig)
        gui.myErrordlg(FigureHandle, ['The number of eigenvalues must be a whole ' ...
            'number of at least 2, or blank.']);
        return;
    end
end

fw = gui.myWaitbar(FigureHandle);
try
    gui.myWaitbar(FigureHandle, fw, false, '', 'Computing the cell-cell covariance spectrum...');
    [T, info] = sc_powerseek(sce.X(:, cellIdx), sce.g, NumGenes=numGenes, ...
        NumEigenvalues=numEig, ...
        ProgressFcn=@(f) gui.myWaitbar(FigureHandle, fw, false, '', ...
        'Scoring genes...', f));
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message);
    return;
end
gui.myWaitbar(FigureHandle, fw);

in_showresult(T, info, cellsName, FigureHandle);
end

function [cellIdx, cellsName] = in_selectcells(src, sce, parentfig)
% Brushed cells win; otherwise one cell type, since differences between
% cell types would dominate the spectrum, or all cells.
cellIdx = [];
cellsName = "";
if isa(src, 'matlab.apps.AppBase') && isprop(src, 'h') && pkg.i_isvalid(src.h)
    brushed = logical(src.h.BrushData(:));
    if any(brushed)
        cellIdx = find(brushed);
        cellsName = sprintf("%d brushed cells", numel(cellIdx));
        return;
    end
end

allCells = sprintf("All cells (n = %d)", sce.NumCells);
types = string(sce.c_cell_type_tx(:));
[typeNames, ~, typeId] = unique(types);
if numel(typeNames) < 2
    cellIdx = (1:sce.NumCells)';
    cellsName = allCells;
    return;
end
counts = accumarray(typeId, 1);
options = [allCells; compose("%s (n = %d)", typeNames, counts)];
[~, order] = sort(counts, "descend");
[indx, tf] = gui.myListdlg(parentfig, options, 'Power-Seek', ...
    order(1) + 1, false, true, [], ...
    ['Run within one cell type: differences between cell types ' ...
    'otherwise dominate the spectrum.']);
if ~tf || isempty(indx), return; end
if indx == 1
    cellIdx = (1:sce.NumCells)';
    cellsName = allCells;
else
    cellIdx = find(typeId == indx - 1);
    cellsName = options(indx);
end
end

function in_showresult(T, info, cellsName, parentfig)
hx = gui.myFigure(parentfig, true);
hFig = hx.FigHandle;
hFig.Name = 'Genes Driving Cell-State Structure (Power-Seek)';
hFig.Position(3) = max(hFig.Position(3), 1200);
hx.addCustomButton('off', @in_callback_showtable, 'list.gif', ...
    'Show the ranked gene table');
hx.addCustomButton('off', @in_callback_savetable, ...
    'floppy-disk-arrow-in.jpg', 'Export the ranked gene table...');

tl = tiledlayout(hFig, 1, 3, TileSpacing="compact", Padding="compact");
title(tl, sprintf('%s: %d cells, %d genes', cellsName, info.NumCells, ...
    info.NumGenes), Interpreter="none");
if isequal(info.HasPowerLaw, false)
    % Without a power-law regime there is no shared-state signal to
    % attribute to genes; the paper's shuffled data rank genes at chance.
    subtitle(tl, ['No clear power-law regime at the top of the spectrum ' ...
        '(log-log fit R^2 < 0.85, or slope under 1.5x the shuffled one). ' ...
        'The ranking is not reliable for these cells.'], Color=[0.8 0 0]);
end
n = info.NumEigenvalues;

% Spectrum against its shuffled null, with the power-law fit.
ax1 = nexttile(tl);
ev = info.Eigenvalues;
loglog(ax1, 1:numel(ev), ev, '.', MarkerSize=8);
hold(ax1, "on");
legendText = {sprintf('Observed (slope %.2f)', info.Exponent)};
if ~isempty(info.EigenvaluesShuffled)
    evs = info.EigenvaluesShuffled;
    loglog(ax1, 1:numel(evs), evs, '.', MarkerSize=8, Color=[0.6 0.6 0.6]);
    legendText{end + 1} = sprintf('Shuffled (slope %.2f)', info.ExponentShuffled);
end
p = polyfit(log((1:n)'), log(ev(1:n)), 1);
loglog(ax1, 1:n, exp(polyval(p, log(1:n))), 'r-', LineWidth=1.5);
legendText{end + 1} = sprintf('Power-law fit, top %d', n);
xline(ax1, n, ':');
hold(ax1, "off");
legend(ax1, legendText, Location="southwest");
xlabel(ax1, 'Rank');
ylabel(ax1, 'Eigenvalue of cell-cell covariance');
title(ax1, 'Eigenvalue spectrum');
box(ax1, "on");

% How the size of the power-law regime was chosen.
ax2 = nexttile(tl);
if isempty(info.RSquared)
    text(ax2, 0.5, 0.5, sprintf('Set by hand: n = %d', n), ...
        HorizontalAlignment="center", Units="normalized");
    axis(ax2, "off");
else
    plot(ax2, info.NumEigenvaluesTried, info.RSquared, '.-');
    xline(ax2, n, 'r:', sprintf('n = %d', n));
    xlabel(ax2, 'Number of top eigenvalues fitted');
    ylabel(ax2, 'R^2 of log-log fit');
    title(ax2, 'Power-law regime (first R^2 peak)');
    box(ax2, "on");
end

% Ranked scores, top genes named.
ax3 = nexttile(tl);
plot(ax3, T.Rank, T.Delta, '.');
yline(ax3, 0, ':');
numLabels = min(10, height(T));
text(ax3, T.Rank(1:numLabels) + 0.01*height(T), T.Delta(1:numLabels), ...
    T.Gene(1:numLabels), FontSize=8, Interpreter="none");
xlabel(ax3, 'Gene rank');
ylabel(ax3, '\Delta (drop in top eigenvalues without the gene)');
title(ax3, sprintf('Ranked genes (%d with \\Delta > 0)', nnz(T.Delta > 0)));
box(ax3, "on");

hx.show(parentfig);

    function in_callback_showtable(~, ~)
        gui.i_openinstead(hFig, @() gui.TableViewerApp(T, hFig, "PowerSeekGenes"));
    end

    function in_callback_savetable(~, ~)
        gui.i_exporttable(T, true, 'Tpowerseek', 'PowerSeekGenes', ...
            [], [], hFig);
    end
end
