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
        sigtags=[sigtags;"-------------------------";
                 "Select scores to make a customized score set..."];
        [indx1,tf1] = listdlg('PromptString','Select a predefined score set',...
             'SelectionMode','single','ListString',...
             sigtags,'ListSize',[250,300]);
        if tf1~=1, return; end
        if contains(sigtags(indx1),'----'), return; end
            
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
                if n~=k
                    gui.gui_waitbar_adv(fw,k/n,listitems{indx2(k)});
                else
                    gui.gui_waitbar_adv(fw,(k-1)/n,listitems{indx2(k)});
                end
                [y]=pkg.e_cellscores(sce.X,sce.g, ...
                    listitems{indx2(k)},1,false);
                Y(:,k)=y(:);
            end
            gui.gui_waitbar_adv(fw);
            T=array2table(Y,'VariableNames', ...
                listitems(indx2),'RowNames', ...
            matlab.lang.makeUniqueStrings(sce.c_cell_id));
            needwait=true;
            gui.i_exporttable(T,needwait);
            %assignin('base','Y',Y);
            %assignin('base','listitems',listitems(indx2));
            %assignin('base','labelx',listitems(indx2));
            
            labelx=listitems(indx2)';
            %gui.gui_waitbar(fw);
            % T=table(Y,'VariableNames', ...
            %     matlab.lang.makeValidName(listitems(indx2)));

    answer = questdlg('Compare between different cell groups?','');
    if ~strcmp(answer,'Yes'), return; end
     allowunique=false;
     [thisc]=gui.i_select1class(sce,allowunique);
     if isempty(thisc), return; end
     [c,cL]=grp2idx(thisc);     
     % assignin('base','thisc',thisc);

     if n==1
         gui.i_violinplot(Y,thisc,labelx,true,[],[]);
         xlabel('Cell group');
         ylabel('Cellular score');
     elseif n==2
         gui.i_violinplot(Y(:,1),thisc,labelx{1},true,[],[]);
         xlabel('Cell group'); ylabel('Cellular score');
         gui.i_violinplot(Y(:,2),thisc,labelx{2},true,[],[]);
         xlabel('Cell group'); ylabel('Cellular score');
    
     elseif n>=3
         P=grpstats(Y,c,'mean');
         assignin('base','P',P);
         assignin('base','labelx',labelx);
         if ~isempty(strfind(labelx{1},')'))
             titlex=extractBefore(labelx{1},strfind(labelx{1},')')+1);
             labelx=extractAfter(labelx,strfind(labelx{1},')')+1);
         else
             titlex='';
         end

         axes_limits=[zeros(1,n); repmat(max(P(:)),1,n)];
         figure;
         spider_plot_R2019b(P,'AxesLabels',labelx, ...
             'AxesPrecision',2,'AxesLimits',axes_limits);
         cL=strrep(cL,'_','\_');
         legend(cL);
         if ~isempty(titlex)
             title(titlex);
         end
     end
end