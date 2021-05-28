function callback_SelectCellsByMarker(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    gsorted=sort(sce.g);
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
        end
        fw=gui.gui_waitbar;
        scex=selectcells(sce,idx);        
        sc_scatter_sce(scex);
        view(ax,bx);
        gui.gui_waitbar(fw);
    end
end