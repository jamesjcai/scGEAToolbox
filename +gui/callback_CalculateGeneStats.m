function callback_CalculateGeneStats(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    fw=gui.gui_waitbar;
    [T]=sc_genestats(sce.X,sce.g)
    gui.gui_waitbar(fw);
    gui.i_exporttable(T);    
end