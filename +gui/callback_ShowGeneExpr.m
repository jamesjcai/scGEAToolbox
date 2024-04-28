function callback_ShowGeneExpr(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    fw=gui.gui_waitbar;
    gui.sc_uitabgrpfig_feaplot(sce, glist, FigureHandle, [axx, bxx]);
    gui.gui_waitbar(fw);
end