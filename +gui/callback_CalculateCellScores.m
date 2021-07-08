function callback_CalculateCellScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    species=questdlg('Use predefined score type or select genes?',...
        'Select Score Type','Select Predefined Score',...
        'Customize Score','Select Predefined Score');
switch species
    case 'Select Predefined Score'
        [~,T]=pkg.e_cellscores(sce.X,sce.g,0);

        listitems=T.ScoreType;
        [indx,tf] = listdlg('PromptString',...
            {'Select Class','',''},...
             'SelectionMode','single','ListString',listitems,'ListSize',[220,300]);
        if ~tf==1, return; end
        [cs]=pkg.e_cellscores(sce.X,sce.g,indx);
        ttxt=T.ScoreType(indx);

    case 'Customize Score'
        ttxt='Customized Score';
        % Pd1=pdcd1 tim3=HAVCR2, tcf1=HNF1A  https://www.nature.com/articles/s41577-019-0221-9    
        % posgcandidates=["PDCD1","HNF1A","HAVCR2","KLRG1","CD44","LY6C","CTLA","ICOS","LAG3"];
        posgcandidates=sce.g(randi(length(sce.g),10,1));
        [posg]=gui.i_selectngenes(sce.g,posgcandidates);
        if isempty(posg)
            helpdlg('No feature genes selected.')
            return;
        end
        posg
        cs=sc_cellscore(sce.X,sce.g,posg);
    otherwise
        return;
end        
    
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Score Value')
    title(strrep(ttxt,'_','\_'))
%     labels = {'Save score values to variable named:'}; 
%     vars = {'CellScores'};
%     values = {cs};
%     export2wsdlg(labels,vars,values);
end


