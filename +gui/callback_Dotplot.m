function callback_Dotplot(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL"). Downstream treats thisc as a
% per-cell label vector, so the composite needs no special handling.
[thisc, clabel] = gui.i_selectnclass(sce,[],[],[],FigureHandle);
if isempty(thisc), return; end
thisc = string(thisc);

% Same group chooser as gui.callback_Violinplot - see gui.i_selectgroupsubset.
[picked, levels] = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
if isempty(picked), return; end
thisc = thisc(picked);

% LEVELS is the order the chooser listed them in, so declining a manual
% order keeps that order rather than re-sorting.
[c, cL, noanswer] = gui.i_reordergroups(thisc, levels, FigureHandle);
if noanswer, return; end

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
    return;
end
[Xt] = gui.i_transformx(sce.X, true, "libsize_log1p", FigureHandle);
if isempty(Xt), return; end
% Subset after normalising, so each cell's values match an all-cells plot.
if ~all(picked), Xt = Xt(:, picked); end
% Not reversed: GUI.I_DOTPLOT draws its first gene on top, so the genes
% read top to bottom in the order picked, as from the heatmap's dot plot
% button. Reversing the list here put the last-picked gene on top.

try
    gui.i_dotplot(Xt, sce.g, c, cL, glist, true, 'Dotplot', FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end
end
