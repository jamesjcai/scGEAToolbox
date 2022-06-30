function callback_TCellExhaustionScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    fw=gui.gui_waitbar;
    try
        cs=pkg.e_cellscores(sce.X,sce.g,"T Cell Exhaustion");
    catch ME
        gui.gui_waitbar(fw,true);
        errordlg(ME.message);
        return;
    end        
    gui.gui_waitbar(fw);
    if ~isempty(cs)
        figure;
        gui.i_stemscatter(sce.s,cs);
        zlabel('Score Value')
        title('T Cell Exhaustion Score')
        gui.i_exporttable(cs,false,'TCellExhaustionScores');
    end
end