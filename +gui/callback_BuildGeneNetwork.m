    function callback_BuildGeneNetwork(src, ~)
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    
    [y, i] = ismember(upper(glist), upper(sce.g));
    if ~all(y), error('xxx'); end
    fprintf("%s\n", glist)
    
    [Xt] = gui.i_transformx(sce.X, true, 5);
    if isempty(Xt), return; end
    
    fw = gui.gui_waitbar;
    x = Xt(i, :);
    A = sc_pcnet(x);
    gui.gui_waitbar(fw);
    sc_grnview(A, glist, [], FigureHandle);
end
