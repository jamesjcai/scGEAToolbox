function callback_CompareGeneBtwCls(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc]=gui.i_select1class(sce);
    if isempty(thisc)
        % warndlg('Grouping variable undefined.');
        return;
    end
    
a={'Library Size','Gene Expression','Predefined Cell Score'};
[indx1,tf1]=listdlg('PromptString',...
    'Select a metric for comparison.',...
    'SelectionMode','single','ListString',a);
if tf1~=1, return; end

switch a{indx1}
    case 'Library Size'
        y=sum(sce.X);
        ttxt='Library Size';
    case 'Gene Expression'
        gsorted=sort(sce.g);
        [indx,tf] = listdlg('PromptString',{'Select a gene',...
                    '',''},'SelectionMode','single',...
                    'ListString',gsorted);
        if tf~=1, return; end
        idx=sce.g==gsorted(indx);
        try
            [Xt]=gui.i_transformx(sce.X);
        catch ME
            errordlg(ME.message);
            return;
        end
        y=full(Xt(idx,:));
        ttxt=sce.g(idx);
    case 'Predefined Cell Score'
        [~,T]=pkg.e_cellscores(sce.X,sce.g,0);
        listitems=T.ScoreType;
        [indx2,tf2] = listdlg('PromptString',...
            {'Select Class','',''},...
             'SelectionMode','single','ListString',...
             listitems,'ListSize',[220,300]);
        if tf2~=1, return; end
        fw=gui.gui_waitbar;
        [y]=pkg.e_cellscores(sce.X,sce.g,indx2);
        ttxt=T.ScoreType(indx2);
        gui.gui_waitbar(fw);        
end

        f = figure('visible','off');
        pkg.i_violinplot(y,thisc);
        title(strrep(ttxt,'_','\_'));
        ylabel('Expression Level');
        movegui(f,'center');
        set(f,'visible','on');    
end