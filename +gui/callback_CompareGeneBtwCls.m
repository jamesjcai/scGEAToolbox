function callback_CompareGeneBtwCls(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [thisc]=gui.i_select1class(sce);
    if isempty(thisc)
        % warndlg('Grouping variable undefined.');
        return;
    end
    gsorted=sort(sce.g);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1
        idx=sce.g==gsorted(indx);
        [Xt]=gui.i_transformx(sce.X);      
        f = figure('visible','off');
        y=full(Xt(idx,:));
        pkg.i_violinplot(y,thisc);
        title(sce.g(idx));
        ylabel('Expression Level')
        movegui(f,'center');
        set(f,'visible','on');
    end
end