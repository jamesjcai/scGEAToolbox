function callback_BuildGeneNetwork(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    
    [glist]=gui.i_selectngenes(sce);
    if isempty(glist), return; end

    [y,i]=ismember(glist,sce.g);
    if ~all(y), error('xxx'); end    
    fprintf("%s\n",glist)
    fw=gui.gui_waitbar;    
    x=sce.X(i,:);
    x=sc_transform(x);
    A=sc_pcnet(x);
    gui.gui_waitbar(fw);            
    sc_grnview(A,glist);
end
