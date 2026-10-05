function callback_Heatmap(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);
% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL"). Downstream treats thisc as a
% per-cell label vector, so the composite needs no special handling.
[thisc, clabel] = gui.i_selectnclass(sce,[],[],[],FigureHandle);
if isempty(thisc), return; end
thisc = string(thisc);

% Same group chooser as gui.callback_Dotplot - see gui.i_selectgroupsubset.
[picked, levels] = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
if isempty(picked), return; end
thisc = thisc(picked);

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
    return;
end
gui.i_heatmap(sce, glist, thisc, FigureHandle, [], levels, picked);
end
