function callback_MacrophagePolarizationIndex(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    helpdlg("Function is under development.");
    return;

    scorename="Macrophage_Polarization_Index";
    cs=pkg.e_cellscores(sce.X,sce.g,scorename);
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Score Value')
    title(scorename)
    
    labels = {'Save score values to variable named:'}; 
    vars = {scorename};
    values = {cs};
    export2wsdlg(labels,vars,values);
end