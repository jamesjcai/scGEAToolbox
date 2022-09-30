function callback_SelectCellsByClass(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    thisc=gui.i_select1class(sce);
    if isempty(thisc), return; end
    [ci,cLi]=grp2idx(thisc);

%     [cLi,idx]=sort(cLi);
%     ci=idx(ci);

    [indxx,tfx] = listdlg('PromptString',{'Select groups'},...
        'SelectionMode','multiple','ListString',string(cLi));
    if tfx==1
        [ax,bx]=view();
        idx=ismember(ci,indxx);

        answer = questdlg('Select or unselect?','','Select', 'Unselect',...
                          'Cancel', 'Select');
        if strcmp(answer, 'Select')
              % do nothing
        elseif strcmp(answer, 'Unselect')            
            idx=~idx;
        else
            return;
        end        
        fw=gui.gui_waitbar;

        scex=selectcells(sce,idx);
        % scex.c=cLi(ci(idx));
        scex.c=ci(idx);
        sc_scatter_sce(scex);
        view(ax,bx);
        gui.gui_waitbar(fw);
    end
    
end