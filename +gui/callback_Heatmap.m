function callback_Heatmap(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
[thisc, ~] = gui.i_select1class(sce,[],[],[],FigureHandle););
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