function callback_SelectCellsByClass(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    thisc=gui.i_select1class(sce);
    if isempty(thisc), return; end
    [ci,cLi]=grp2idx(thisc);

    [indxx,tfx] = listdlg('PromptString',{'Select groups',...
    '',''},'SelectionMode','multiple','ListString',string(cLi));
    if tfx==1
        [ax,bx]=view();
        fw=gui.gui_waitbar;
        idx=ismember(ci,indxx);
        scex=selectcells(sce,idx);
        scex.c=cLi(ci(idx));
        sc_scatter_sce(scex);
        view(ax,bx);
        gui.gui_waitbar(fw);
    end
    
end