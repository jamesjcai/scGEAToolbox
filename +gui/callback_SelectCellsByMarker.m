function callback_SelectCellsByMarker(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    
    % gsorted=sort(sce.g);

answer1 = questdlg('Use single or mulitple markers?',...
    'Single/Multiple Markers','Single','Multiple','Cancel','Single');

switch answer1
    case 'Single'
        do_single;
    case 'Multiple'
        do_multiple;
    otherwise
        return;
end

function do_multiple
    [glist]=gui.i_selectngenes(sce);
    if ~isempty(glist)
       [y,i]=ismember(upper(glist),upper(sce.g));
       if ~all(y), error('Unspecific running error.'); end
        ix=sum(sce.X(i,:)>0,1)==length(i);
        if ~any(ix)
            helpdlg('No cells expressing all selected markers.','');
            return;
        end
       
       answer = questdlg('Extract markers+ or markers- cells?',...
           'Positive or Negative',...
            'Markers+','Markers-','Cancel','Markers+');
       switch answer
           case 'Markers+'
                idx=ix;
           case 'Markers-'
                idx=~ix;
           case 'Cancel'
               return;
           otherwise
               return;
       end
        [ax,bx]=view();
        fw=gui.gui_waitbar;
        scex=selectcells(sce,idx);        
        sc_scatter_sce(scex);
        view(ax,bx);
        gui.gui_waitbar(fw);

    end
end

function do_single
    [gsorted]=gui.i_sortgenenames(sce);
    if isempty(gsorted), return; end

    [indx,tf] = listdlg('PromptString',{'Select a gene','',''},...
                'SelectionMode','single','ListString',gsorted);
    if tf==1
        [ax,bx]=view();
        tg=gsorted(indx);
        c = sce.X(sce.g == tg, :);
        answer = questdlg(sprintf('Extract %s+ or %s- cells?',tg,tg),'Positive or Negative',...
            sprintf('%s+',tg),sprintf('%s-',tg),'Cancel',sprintf('%s+',tg));
        if strcmp(answer,sprintf('%s+',tg))
            idx=c>0;
        elseif strcmp(answer,sprintf('%s-',tg))
            idx=c==0;
        elseif strcmp(answer,'Cancel')
           return;
        else
            return;
        end
        fw=gui.gui_waitbar;
        scex=selectcells(sce,idx);        
        sc_scatter_sce(scex);
        view(ax,bx);
        gui.gui_waitbar(fw);
    end
end

end