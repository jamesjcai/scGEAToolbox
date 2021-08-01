function callback_CalculateCellScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    actiontype=questdlg('Select a predefined score or define a new score?',...
        'Select/Define Score','Select Predefined Score',...
        'Define New Score','Select Predefined Score');
switch actiontype
    case 'Select Predefined Score'
        [~,T]=pkg.e_cellscores(sce.X,sce.g,0);

        listitems=T.ScoreType;
        [indx,tf] = listdlg('PromptString',...
            {'Select Class','',''},...
             'SelectionMode','single','ListString',listitems,'ListSize',[220,300]);
        if ~tf==1, return; end
        fw=gui.gui_waitbar;
        [cs]=pkg.e_cellscores(sce.X,sce.g,indx);
        ttxt=T.ScoreType(indx);
        gui.gui_waitbar(fw);
        
    case 'Define New Score'
        ttxt='Customized Score';
        % Pd1=pdcd1 tim3=HAVCR2, tcf1=HNF1A  https://www.nature.com/articles/s41577-019-0221-9    
        % posgcandidates=["PDCD1","HNF1A","HAVCR2","KLRG1","CD44","LY6C","CTLA","ICOS","LAG3"];
        %posgcandidates=sce.g(randi(length(sce.g),10,1));
        [posg]=gui.i_selectngenes(sce.g);
        if isempty(posg)
            helpdlg('No feature genes selected.')
            return;
        end
        fw=gui.gui_waitbar;
        ttxt=sprintf('Customized Score\n%s',sprintf('%s+',posg));
        
        cs=sc_cellscore(sce.X,sce.g,posg);
        gui.gui_waitbar(fw);
    otherwise
        return;
end        
    
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Score Value')
    title(strrep(ttxt,'_','\_'))
    
    
    answer2 = questdlg(sprintf('Score has been computed.\nCompare the score across cell classes?'),'Continue?');
    switch answer2
        case 'Yes'
    
        otherwise
            return;
    end
    [thisc,clabel]=gui.i_select1class(sce);
    if isempty(thisc)   % || numel(unique(thisc))==1
        errordlg('Undefined');
        return;
    end
    figure;
    pkg.i_violinplot_groupordered(cs,thisc);
    ylabel(strrep(ttxt,'_','\_'))
    xlabel(clabel);
    
%         f = figure('visible','off');
%         y=full(Xt(idx,:));
%         pkg.i_violinplot(y,thisc);
%         title(sce.g(idx));
%         ylabel('Expression Level')
%         movegui(f,'center');
%         set(f,'visible','on');
    
%     labels = {'Save score values to variable named:'}; 
%     vars = {'CellScores'};
%     values = {cs};
%     export2wsdlg(labels,vars,values);
end


