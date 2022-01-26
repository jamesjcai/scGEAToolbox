function callback_TCellExhaustionScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    fw=gui.gui_waitbar;
    cs=pkg.e_cellscores(sce.X,sce.g,"T_Cell_Exhaustion");
    gui.gui_waitbar(fw);
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Score Value')
    title('T Cell Exhaustion Score')
    gui.i_exporttable(cs,false,'TCellExhaustionScores');

%     if ~(ismcc || isdeployed)
%         labels = {'Save score values to variable named:'}; 
%         vars = {'TCellExhaustionScores'};
%         values = {cs};
%         export2wsdlg(labels,vars,values);
%     else
%         gui.i_exporttable(cs,false,'TCellExhaustionScores');
%     end
end