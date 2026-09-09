function [needupdate] = callback_InferCNV(src, ~)
%CALLBACK_INFERCNV Infer copy number variation and call malignant cells.
%
%   Runs SC_INFERCNV followed by SC_MALIGNSCORE on the current SCE. The user
%   nominates a group of cells that are assumed karyotypically normal; every
%   result is relative to them, so the choice matters more than any other
%   setting here.
%
%   Adds two cell attributes, viewable under View -> Cell State (Ctrl + T):
%     malignancy_score  copy-number deviation from the reference cells
%     malignancy_type   "malignant" / "nonMalignant", when the score splits
%
%   See also SC_INFERCNV, SC_MALIGNSCORE, GUI.CALLBACK_CELLCYCLEPOTENCY.

needupdate = false;
[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~gui.gui_showrefinfo('InferCNV', FigureHandle), return; end

% ---- Genome ----
speciestag = gui.i_selectspecies(2, false, FigureHandle);
if isempty(speciestag), return; end
if strcmpi(speciestag, 'mouse')
    genome = 'mm10';
else
    answer = gui.myQuestdlg(FigureHandle, 'Which genome build?', ...
        'Select Genome', {'hg38', 'hg19'}, 'hg38');
    if isempty(answer), return; end
    genome = answer;
end

% ---- Reference cells ----
[thisc, clabel] = gui.i_select1class(sce, false, ...
    'Which grouping variable identifies the normal cells?', ...
    'Cell Type', FigureHandle);
if isempty(thisc), return; end

[isref, levels] = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
if isempty(isref), return; end

isref = logical(isref(:)).';
nref = sum(isref);
if nref == 0
    gui.myWarndlg(FigureHandle, ...
        'No reference cells were selected. Pick at least one group.');
    return;
end
if all(isref)
    gui.myWarndlg(FigureHandle, ...
        ['Every cell was selected as reference, leaving nothing to test ' ...
        'against. Pick only the groups you believe are normal.']);
    return;
end
if nref < 50
    answer = gui.myQuestdlg(FigureHandle, sprintf( ...
        ['Only %d reference cells were selected. The reference profile is ' ...
        'the baseline for every other cell, and a noisy one hides real ' ...
        'events. A few hundred is preferable. Continue anyway?'], nref));
    if ~strcmp(answer, 'Yes'), return; end
end

% ---- Run ----
fw = gui.myWaitbar(FigureHandle);
try
    [cnv, T] = sc_infercnv(sce.X, sce.g, isref, Genome=genome, Verbose=false);

    % Pooling a cell with its neighbours sharpens the score, because dropout
    % makes single-cell CNV profiles noisy while cells of one clone share a
    % karyotype.
    %
    % Test for a real embedding with I_CHECKEXISTINGEMBED, not with
    % ~isempty(sce.s): a freshly constructed SingleCellExperiment already
    % carries random coordinates in .s, so the obvious check passes on data
    % that has never been embedded and would pool each cell with a random
    % handful of others -- quietly making the score worse rather than
    % better. STRUCT_CELL_EMBEDDINGS is only filled by an actual embedding.
    A = [];
    if ~isempty(gui.i_checkexistingembed(sce)) && sce.NumCells <= 50000
        A = sc_knngraph(sce.s(:, 1:min(3, size(sce.s, 2))), 10, false);
    end
    if isempty(A)
        [score, label, info] = sc_malignscore(cnv, isref, Verbose=false);
    else
        [score, label, info] = sc_malignscore(cnv, isref, Adjacency=A, ...
            Verbose=false);
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'InferCNV');
    return;
end
gui.myWaitbar(FigureHandle, fw);

% ---- Store ----
sce.setCellAttribute('malignancy_score', score);
sce.setCellAttribute('malignancy_type', string(label));
gui.myGuidata(FigureHandle, sce, src);
needupdate = true;

nmal = sum(label == "malignant");
if info.hasSplit
    verdict = sprintf('%d of %d cells called malignant.', nmal, sce.NumCells);
else
    verdict = ['The scores are unimodal, so no malignant population was ' ...
        'separated and every cell is labelled nonMalignant. That is a ' ...
        'finding, not a failure: it is what a sample with no aneuploid ' ...
        'cells looks like.'];
end
gui.myHelpdlg(FigureHandle, sprintf( ...
    ['%s\n\nReference: %d cells from %s (%s).\n\n' ...
    'malignancy_score and malignancy_type added. To see them, use ' ...
    'View -> Cell State (Ctrl + T).'], ...
    verdict, nref, clabel, strjoin(cellstr(levels(:).'), ', ')), 'InferCNV');

% ---- Optional heatmap ----
answer = gui.myQuestdlg(FigureHandle, ...
    'Show the CNV heatmap (genes along the genome by cells)?');
if strcmp(answer, 'Yes')
    in_cnvheatmap(cnv, T, label, isref, FigureHandle);
end

end


function in_cnvheatmap(cnv, T, label, isref, parentfig)
% Genome on the x axis, cells on the y axis, grouped reference first so the
% contrast between the baseline and anything aneuploid is visible at a
% glance. Rows are thinned rather than plotted in full: a few thousand
% screen pixels cannot show more.

maxCells = 3000;
order = [find(isref(:).' & label(:).' ~= "malignant"), ...
    find(~isref(:).' & label(:).' ~= "malignant"), ...
    find(label(:).' == "malignant")];
if numel(order) > maxCells
    order = order(round(linspace(1, numel(order), maxCells)));
end

M = double(cnv(:, order)).';
nRef = sum(isref(order));
nNon = sum(label(order) ~= "malignant");

hFig = figure('Name', 'InferCNV', 'Color', 'w', 'NumberTitle', 'off');
if ~isempty(parentfig)
    gui.i_movegui2parent(hFig, parentfig);
end
ax = axes(hFig);

imagesc(ax, M);
clim(ax, [0.85 1.15]);
colormap(ax, i_cnvcolormap());
cb = colorbar(ax);
cb.Label.String = 'relative copy number';

% Chromosome boundaries and centred labels.
edges = [1; find(diff(T.Chr) ~= 0) + 1; height(T) + 1];
hold(ax, 'on');
for k = 2:numel(edges)-1
    xline(ax, edges(k) - 0.5, 'k-', 'Alpha', 0.35, 'LineWidth', 0.5);
end
% Separate the reference block, then the rest of the non-malignant cells.
yline(ax, nRef + 0.5, 'b-', 'LineWidth', 1.2, 'Alpha', 0.8);
if nNon > nRef
    yline(ax, nNon + 0.5, 'r-', 'LineWidth', 1.2, 'Alpha', 0.8);
end
hold(ax, 'off');

chrs = unique(T.Chr, 'stable');
ax.XTick = (edges(1:end-1) + edges(2:end) - 1) / 2;
ax.XTickLabel = string(chrs);
ax.XAxis.FontSize = 7;
xlabel(ax, 'chromosome');
ylabel(ax, sprintf('cells (n=%d shown; blue = end of reference)', numel(order)));
title(ax, 'Inferred copy number');
end


function cmap = i_cnvcolormap()
% Blue for loss, white for neutral, red for gain, with white at exactly the
% midpoint so that no-change reads as no-change.
n = 128;
half = linspace(0, 1, n).';
cmap = [ [half, half, ones(n, 1)]; [ones(n, 1), flipud(half), flipud(half)] ];
end
