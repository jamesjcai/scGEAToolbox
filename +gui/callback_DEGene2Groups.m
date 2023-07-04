function callback_DEGene2Groups(src,~)
    
    isatac=false;

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    

    [i1,i2,cL1,cL2]=gui.i_select2grps(sce);
    if length(i1)==1 || length(i2)==1, return; end

    answer = questdlg('Which method?',...
        'Select Method','Wilcoxon rank-sum test 🐇',...
        'DESeq 2 (R required) 🐢',...
        'MAST (R required) 🐢', ...
        'Wilcoxon rank-sum test 🐇');
    
    if strcmpi(answer,'Wilcoxon rank-sum test 🐇')
        methodtag="ranksum";
    elseif strcmpi(answer,'DESeq 2 (R required) 🐢')
        methodtag="deseq2";
%         if ~(ismcc || isdeployed)
%             if ~exist('nbintest.m', 'file')
%                 errordlg('This option requires Bioinformatics toolbox.');
%                 return;
%             end
%         end
    elseif strcmpi(answer,'MAST (R required) 🐢')
        methodtag="mast";
        if isempty(pkg.FindRpath)
            warndlg('This function requires R environment.')
            return;
        end
    else
        return;
    end
    try
        switch methodtag
            case 'ranksum'
                %fw=gui.gui_waitbar;
                T=sc_deg(sce.X(:,i1),sce.X(:,i2), ...
                    sce.g,1,true);
            case 'deseq2'
                [ok]=gui.i_confirmscript('DE analysis (DESeq2)', ...
                    'R_DESeq2','r');
                if ~ok, return; end
                fw=gui.gui_waitbar;
                T=run.r_DESeq2(sce.X(:,i1),sce.X(:,i2),sce.g);
                gui.gui_waitbar(fw);
            case 'mast'
                [ok]=gui.i_confirmscript('DE analysis (MAST)', ...
                    'R_MAST','r');
                if ~ok, return; end
                fw=gui.gui_waitbar;
                T=run.r_MAST(sce.X(:,i1),sce.X(:,i2),sce.g);
                gui.gui_waitbar(fw);
        end

    catch ME
        gui.gui_waitbar(fw,true);
        errordlg(ME.message);
        return;
    end


% figure;
% gui.i_volcanoplot(T);
% title(sprintf('%s vs. %s', ...
%     matlab.lang.makeValidName(string(cL1)),matlab.lang.makeValidName(string(cL2))));


% T2=T;
% T2.avg_log2FC(T.avg_log2FC>10)=10;
% T2.avg_log2FC(T.avg_log2FC<-10)=-10;
% T2.p_val_adj(T.p_val_adj<1e-50)=1e-50;
% idx=(T2.avg_log2FC>1 | T2.avg_log2FC<-1) & -log10(T2.p_val_adj)>2;
% scatter(T2.avg_log2FC,-log10(T2.p_val_adj),10,idx+1);
% xline(0); xline(-1); xline(1);
% yline(2);
% colormap(gca,lines(2));

    % mavolcanoplot(sce.X(:,i1),sce.X(:,i2),T.p_val_adj,'Labels',T.gene)

    try
        T = sortrows(T,'p_val_adj','ascend');
        T = sortrows(T,'pct_1','ascend');
        T = sortrows(T,'pct_2','descend');
        T = sortrows(T,'avg_log2FC','ascend');
    catch ME
        warning(ME.message);
    end

try

% avg_1 = mean(X,2);
% avg_2 = mean(Y,2);
% pct_1 = sum(X>0,2)./size(X,2);
% pct_2 = sum(Y>0,2)./size(Y,2);

if contains(T.Properties.VariableNames{5},'avg_1')
    T.Properties.VariableNames{5}=sprintf('%s_%s', ...
        T.Properties.VariableNames{5}, ...
        matlab.lang.makeValidName(string(cL1)));
end

if contains(T.Properties.VariableNames{6},'avg_2')
    T.Properties.VariableNames{6}=sprintf('%s_%s', ...
        T.Properties.VariableNames{6}, ...
        matlab.lang.makeValidName(string(cL2)));
end

if contains(T.Properties.VariableNames{7},'pct_1')
    T.Properties.VariableNames{7}=sprintf('%s_%s', ...
        T.Properties.VariableNames{7}, ...
        matlab.lang.makeValidName(string(cL1)));
end

if contains(T.Properties.VariableNames{8},'pct_2')
    T.Properties.VariableNames{8}=sprintf('%s_%s', ...
        T.Properties.VariableNames{8}, ...
        matlab.lang.makeValidName(string(cL2)));
end
catch ME
    warning(ME.message);
end

    outfile=sprintf('%s_vs_%s', ...
        matlab.lang.makeValidName(string(cL1)),matlab.lang.makeValidName(string(cL2)));
    if isatac, T.gene="chr"+T.gene; end

    [filetype,filesaved]=gui.i_exporttable(T,true,'T', outfile);
    tf=0;
    if ~(ismcc || isdeployed) && strcmp(filetype,'Workspace')
        [Tup,Tdn]=pkg.e_processDETable(T,true);
        labels = {'Save DE results (selected up-regulated) to variable named:',...
            'Save DE results (selected down-regulated) to variable named:'}; 
        vars = {'Tup','Tdn'}; values = {Tup,Tdn};
        [~,tf]=export2wsdlg(labels,vars,values);
    end
   
    if ~isempty(filesaved)
        if strcmp(filetype,'Excel file')
            answer = questdlg('Save up- and down-regulated genes to seperate sheets?');
            if strcmp(answer,'Yes')            
                [Tup,Tdn]=pkg.e_processDETable(T,true);
                % strcmp(extractAfter(filesaved,strlength(filesaved)-4),'xlsx')
                writetable(Tup,filesaved,"FileType","spreadsheet",'Sheet','Up-regulated');
                writetable(Tdn,filesaved,"FileType","spreadsheet",'Sheet','Down-regulated'); 
                waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
                %writetable(Tup,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet',);
                %writetable(Tdn,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet');
            end
        elseif strcmp(filetype,'Text file')
            % strcmp(extractAfter(filesaved,strlength(filesaved)-3),'txt')
            answer = questdlg('Save up- and down-regulated genes to seperate files?');
            if strcmp(answer,'Yes')
                [Tup,Tdn]=pkg.e_processDETable(T,true);
                [~,~]=gui.i_exporttable(Tup,true,'Tup','Upregulated','Text file');
                [~,~]=gui.i_exporttable(Tdn,true,'Tdn','Downregulated','Text file');
            end
        end
    end



%    if ~(ismcc || isdeployed)
        % answer = questdlg('Save up- and down-regulated genes for enrichment analysis?');
        % 
        % if strcmp(answer,'Yes')
        %     [Tup,Tdn]=pkg.e_processDETable(T,true);
        %     tf=0;
        %     if ~(ismcc || isdeployed) && strcmp(filetype,'Workspace')
        %         labels = {'Save DE results (selected up-regulated) to variable named:',...
        %             'Save DE results (selected down-regulated) to variable named:'}; 
        %         vars = {'Tup','Tdn'}; values = {Tup,Tdn};
        %         [~,tf]=export2wsdlg(labels,vars,values);
        %     end
        % 
        %     if ~isempty(filesaved)
        %         if strcmp(filetype,'Excel file')
        %             % strcmp(extractAfter(filesaved,strlength(filesaved)-4),'xlsx')
        %             writetable(Tup,filesaved,"FileType","spreadsheet",'Sheet','Up-regulated');
        %             writetable(Tdn,filesaved,"FileType","spreadsheet",'Sheet','Down-regulated'); 
        %             %waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
        %             %writetable(Tup,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet',);
        %             %writetable(Tdn,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet');
        %         elseif strcmp(filetype,'Text file')
        %             % strcmp(extractAfter(filesaved,strlength(filesaved)-3),'txt')
        %             [~,filesaved1]=gui.i_exporttable(Tup,true,'Tup');
        %             if ~isempty(filesaved1)
        %                 waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved1),''));
        %             end
        %             [~,filesaved2]=gui.i_exporttable(Tdn,true,'Tdn');
        %             if ~isempty(filesaved2)
        %                 waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved2),''));
        %             end
        %         end
        %     end

            if tf==1
                disp('To run enrichment analysis, type:')
                disp('run.web_Enrichr(Tup.gene(1:200))')
                disp('run.web_Enrichr(Tdn.gene(1:200))')
            end
        %end
        %return;
    
    
        answer = questdlg('Run enrichment analysis with top K (=100 by default) up-regulated DE genes?');
        if strcmp(answer,'Yes')
            gui.i_enrichtest(Tup.gene(1:min(numel(Tup.gene),100)));
        end
        
        answer = questdlg('Run enrichment analysis with top K (=100 by default) down-regulated DE genes?');
        if strcmp(answer,'Yes')
            gui.i_enrichtest(Tdn.gene(1:min(numel(Tdn.gene),100)));
        end

%    end
    
%     pause(3);
%     
%     %Tf=run.fgsea(T.gene);
%     %pkg.e_fgseanet(Tf);
    
%     answer = questdlg('Violin plots (top 16 DE genes)?');
%     if strcmp(answer,'Yes')
%         figure;
%         for k=1:16
%             subplot(4,4,k)
%             i=sce.g==T.gene(k);
%             pkg.i_violinplot(log2(1+sce.X(i,:)),...
%                 sce.c_batch_id);
%             title(T.gene(k));
%             ylabel('log2(UMI+1)')
%         end
%     end

%     answer = questdlg('Volcano plot?');
%     if strcmp(answer,'Yes')
%         figure;
%         gui.i_volcanoplot(T,isok);
%     end
    
end