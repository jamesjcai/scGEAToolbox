function callback_CleanEmbeddingPlot(src, ~)
%CALLBACK_CLEANEMBEDDINGPLOT Copy the main embedding plot without axes, for figures.
%
%   Behind Plots > Clean Embedding Plot for Figures. The main window's plot
%   is copied as shown - colours, group labels - with the axes replaced by
%   two (or three) short labelled arrows.
%
%   Narrowing to a subset of cells is an opt-in first step, through the
%   group chooser the Dotplot, Heatmap and violin plots use (see
%   gui.i_selectgroupsubset). The copy then keeps only those cells and the
%   group labels on them, framed to fit.
%
%   See also GUI.I_BAREDRPLOT.

[FigureHandle, sce] = gui.gui_getfigsce(src);

picked = [];
answer = gui.myQuestdlg(FigureHandle, ...
    "Plot all cells, or only cells in selected groups?", "", ...
    {'All Cells', 'Selected Groups', 'Cancel'}, 'All Cells');
switch answer
    case 'All Cells'
        % Keep every cell.
    case 'Selected Groups'
        [thisc, clabel] = gui.i_selectnclass(sce, false, [], [], FigureHandle);
        if isempty(thisc), return; end
        picked = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
        if isempty(picked), return; end
    otherwise
        return;
end

% Label the axis arrows with the method behind SCE.S; ask only when that
% cannot be told from the stored embeddings.
label = pkg.i_embeddinglabel(sce);
if strlength(label) == 0
    label = gui.myQuestdlg(FigureHandle, 'Select embedding method label.', ...
        '', {'tSNE', 'UMAP', 'PHATE'}, 'tSNE');
    if isempty(label), return; end
end
a = colormap(src.UIAxes);
ax2 = gui.i_baredrplot(src.UIAxes, [], char(label), FigureHandle, picked);
colormap(ax2, a);
end
