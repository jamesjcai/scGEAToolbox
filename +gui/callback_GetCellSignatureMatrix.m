function callback_GetCellSignatureMatrix(src,~)
%     answer = questdlg(['This function ' ...
%         'calculates selected signature scores for each ' ...
%         'cell. You will get a signature matrix for cells.' ...
%         ' Continue?'],'');
%     if ~strcmp(answer,'Yes'), return; end

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    preselected=[];

    [~,T]=pkg.e_cellscores([],[],0);
    sigtags=unique(string(T.SignatureTag));
    sigtags=sigtags(strlength(sigtags)>0);
    if ~isempty(sigtags)
        sigtags=[sigtags;"Select score set..."];
        [indx1,tf1] = listdlg('PromptString','Select a signature',...
             'SelectionMode','single','ListString',...
             sigtags,'ListSize',[250,300]);
        if tf1~=1, return; end
            
        idx=T.SignatureTag==sigtags(indx1);
        if any(idx)
            [~,ix]=natsort(T.ScoreType);
            preselected=idx(ix);
        end
    end
            listitems=natsort(T.ScoreType);
            [indx2,tf2] = listdlg('PromptString','Select Scores',...
                 'SelectionMode','multiple','ListString',...
                 listitems,'ListSize',[320,300], ...
                 'InitialValue',find(preselected));
            if tf2~=1, return; end

            n=length(indx2);
            Y=zeros(sce.NumCells,n);

            fw=gui.gui_waitbar_adv;
            for k=1:n
                gui.gui_waitbar_adv(fw,k/n,listitems{indx2(k)});
                [y]=pkg.e_cellscores(sce.X,sce.g, ...
                    listitems{indx2(k)},1,false);
                Y(:,k)=y(:);
            end
            gui.gui_waitbar_adv(fw);
            T=array2table(Y,'VariableNames', ...
                listitems(indx2),'RowNames', ...
            matlab.lang.makeUniqueStrings(sce.c_cell_id));
            gui.i_exporttable(T);
            %assignin('base','Y',Y);
            %assignin('base','listitems',listitems(indx2));
            %gui.gui_waitbar(fw);
            % T=table(Y,'VariableNames', ...
            %     matlab.lang.makeValidName(listitems(indx2)));
end