function callback_CalculateGeneStats(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [T]=sc_genestats(sce.X,sce.g)
    gui.i_exporttable(T);    
end