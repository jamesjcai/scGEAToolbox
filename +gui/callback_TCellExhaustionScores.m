function callback_TCellExhaustionScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    cs=sc_cellscore(sce.X,sce.g);
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Exhaustion Score')
    title('T Cell Exhaustion Score')
end