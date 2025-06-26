function callback_Heatmap(src, ~)

if isa(src, 'matlab.apps.AppBase')    
    [FigureHandle, sce] = xui.gui_getfigsce(src);
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end
[thisc, ~] = gui.i_select1class(sce,[],[],[],FigureHandle);
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