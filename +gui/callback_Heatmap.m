function callback_Heatmap(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);
% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL"). Downstream treats thisc as a
% per-cell label vector, so the composite needs no special handling.
[thisc, ~] = gui.i_selectnclass(sce,[],[],[],FigureHandle);
if isempty(thisc), return; end

% [c, cL, noanswer] = gui.i_reordergroups(thisc);
% if noanswer, return; end

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
    return;
end
gui.i_heatmap(sce, glist, thisc, FigureHandle);
end
