function [requirerefresh,highlightindex]=callback_AnnotateSubGroup(src,~)

requirerefresh=false;
highlightindex=[];

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    thisc=gui.i_select1class(sce);
    if isempty(thisc), return; end
    [ci,cLi]=grp2idx(thisc);

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

        gui.i_gscatter3(scex.s,scex.c);

        % sc_scatter_sce(scex);
        view(ax,bx);
        gui.gui_waitbar(fw);
        highlightindex=idx;
        requirerefresh=true;
    end
    
end