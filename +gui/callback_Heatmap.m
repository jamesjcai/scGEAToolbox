function callback_Heatmap(src, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
[thisc, ~] = gui.i_select1class(sce);
if isempty(thisc), return; end

% [c, cL, noanswer] = gui.i_reordergroups(thisc);
% if noanswer, return; end

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    helpdlg('No gene selected.', '');
    return;
end
gui.i_heatmap(sce, glist, thisc, FigureHandle);
end