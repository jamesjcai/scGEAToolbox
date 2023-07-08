function callback_GetCellSignatureMatrix(src,~)
%     answer = questdlg(['This function ' ...
%         'calculates selected signature scores for each ' ...
%         'cell. You will get a signature matrix for cells.' ...
%         ' Continue?'],'');
%     if ~strcmp(answer,'Yes'), return; end

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);



            [~,T]=pkg.e_cellscores([],[],0);
            listitems=natsort(T.ScoreType);
            [indx2,tf2] = listdlg('PromptString','Select Scores',...
                 'SelectionMode','multiple','ListString',...
                 listitems,'ListSize',[320,300]);
            if tf2~=1, return; end

            n=length(indx2);
            Y=zeros(sce.NumCells,n);

            %fw=gui.gui_waitbar;
            for k=1:n
                [y]=pkg.e_cellscores(sce.X,sce.g, ...
                    listitems{indx2(k)},1);
                %ttxt=T.ScoreType(indx2);
                Y(:,k)=y(:);
            end
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