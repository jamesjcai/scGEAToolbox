function callback_TCellExhaustionScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    
       
    % Pd1=pdcd1 tim3=HAVCR2, tcf1=HNF1A  https://www.nature.com/articles/s41577-019-0221-9    
    posgcandidates=["PDCD1","HNF1A","HAVCR2","KLRG1","CD44","LY6C","CTLA","ICOS","LAG3"];
    [posg]=gui.i_selectngenes(sce.g,posgcandidates);
    if isempty(posg)
        return;
    end
    %negg=["HAVCR2"];
    cs=sc_cellscore(sce.X,sce.g,posg);
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Exhaustion Score')
    title('T Cell Exhaustion Score')
end