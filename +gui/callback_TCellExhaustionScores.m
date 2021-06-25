function callback_TCellExhaustionScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    cs=pkg.e_cellscores(sce.X,sce.g,"T_Cell_Exhaustion");
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Score Value')
    title('T Cell Exhaustion Score')
            labels = {'Save score values to variable named:'}; 
        vars = {'TCellExhaustionScores'};
        values = {cs};
        export2wsdlg(labels,vars,values);
end