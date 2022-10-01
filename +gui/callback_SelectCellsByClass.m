function callback_SelectCellsByClass(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    
    [ptsSelected]=gui.i_select1classcells(sce,true);
    if isempty(ptsSelected), return; end

%{
    thisc=gui.i_select1class(sce);
    if isempty(thisc), return; end
    [ci,cLi]=grp2idx(thisc);

    % [cLi,idx]=sort(cLi);
    % ci=idx(ci);
    [indxx,tfx] = listdlg('PromptString',{'Select groups'},...
        'SelectionMode','multiple','ListString',string(cLi));
    if tfx==1
        [ax,bx]=view();
        ptsSelected=ismember(ci,indxx);

        answer = questdlg('Select or unselect?','','Select', 'Unselect',...
                          'Cancel', 'Select');
        if strcmp(answer, 'Select')
              % do nothing
        elseif strcmp(answer, 'Unselect')            
            ptsSelected=~ptsSelected;
        else
            return;
        end

%}

[ax,bx]=view();
        fw=gui.gui_waitbar;
        scex=selectcells(sce,ptsSelected);
        % scex.c=cLi(ci(idx));
        scex.c=sce.c(ptsSelected);
        sc_scatter_sce(scex);
        view(ax,bx);
        gui.gui_waitbar(fw);
end
    
