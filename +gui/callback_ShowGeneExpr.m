function callback_ShowGeneExpr(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    fw=gui.gui_waitbar;

    n = length(glist);
    y = cell(n,1);
    for k=1:n
        y{k} = sce.X(sce.g == glist(k), :);
    end

    gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
    gui.gui_waitbar(fw);
end