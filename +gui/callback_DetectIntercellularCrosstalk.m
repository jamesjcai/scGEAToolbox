function callback_DetectIntercellularCrosstalk(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    
    [c,cL]=grp2idx(sce.c_cell_type_tx);    
    [idx]=gui.gui_selmultidlg(cL);    
    if isempty(idx), return; end
    if length(idx)<2
        warndlg('Need at least 2 cell types');
        return;
    end
    selected=ismember(c,idx);    
    fw=gui.gui_waitbar;
    sce=sce.selectcells(selected);
    [Tok,OUT]=run.talklr(sce);
    gui.gui_waitbar(fw);
    labels = {'Save T to variable named:', ...
              'Save OUT to variable named:'};
    vars = {'T', 'OUT'};
    values = {Tok, OUT};
    export2wsdlg(labels, vars, values);
    disp('To get more crosstalk network, run gui.i_crosstalkgraph(OUT,k)')
end
