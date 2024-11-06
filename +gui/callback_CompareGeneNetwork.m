function callback_CompareGeneNetwork(src, ~)
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    
    [i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false);
    if isscalar(i1) || isscalar(i2)
        if i1 == 0 || i2 == 0, return; end
    end
    
    [glist] = gui.i_selectngenes(sce,[],FigureHandle);
    if isempty(glist), return; end
    
    [y, i] = ismember(glist, sce.g);
    if ~all(y), error('Selected gene(s) not in the gene list of data.'); end
    fprintf("%s\n", glist)
    
    [Xt] = gui.i_transformx(sce.X, true, 5);
    if isempty(Xt), return; end
    
    fw=gui.gui_waitbar;
    x1 = Xt(i, i1);
    x2 = Xt(i, i2);
    A1 = sc_pcnet(x1, 3, false, true, false);
    A2 = sc_pcnet(x2, 3, false, true, false);
    gui.gui_waitbar(fw);
    pause(1)
    stitle = sprintf('%s vs. %s', cL1{1}, cL2{1});
    sc_grnview2(A1, A2, glist, stitle, FigureHandle);
end
