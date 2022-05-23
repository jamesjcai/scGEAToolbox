function callback_DEGene2Groups(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    

    [i1,i2]=gui.i_select2grps(sce);
    if length(i1)==1 || length(i2)==1, return; end

    answer = questdlg('Which method?',...
        'Select Method','Wilcoxon rank-sum test 🐇',...
        'MAST (R required) 🐢','Wilcoxon rank-sum test 🐇');
    
    if strcmpi(answer,'Wilcoxon rank-sum test 🐇')
        methodtag="ranksum";
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
                fw=gui.gui_waitbar;
                T=sc_deg(sce.X(:,i1),sce.X(:,i2),sce.g);
            case 'mast'
                [ok]=gui.i_confirmscript('DE analysis (MAST)', ...
                    'R_MAST','r');
                if ~ok, return; end
                fw=gui.gui_waitbar;
                T=run.MAST(sce.X(:,i1),sce.X(:,i2),sce.g);
        end
    catch ME
        gui.gui_waitbar(fw,true);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);

    try
        T = sortrows(T,'p_val_adj','ascend');
        T = sortrows(T,'pct_1','ascend');
        T = sortrows(T,'pct_2','descend');
        T = sortrows(T,'avg_log2FC','ascend');
    catch ME
        warning(ME.message);
    end
    gui.i_exporttable(T,true);

    if ~(ismcc || isdeployed)  
        answer = questdlg('Save up- and down-regulated genes for enrichment analysis?');
        if strcmp(answer,'Yes')
            [Tup,Tdn]=pkg.e_processDETable(T);
            labels = {'Save DE results (selected up-regulated) to variable named:',...
                'Save DE results (selected down-regulated) to variable named:'}; 
            vars = {'Tup','Tdn'}; values = {Tup,Tdn};
            [~,tf]=export2wsdlg(labels,vars,values);
            if tf==1
                disp('To run enrichment analysis, type:')
                disp('run.Enrichr(Tup.gene(1:200))')
                disp('run.Enrichr(Tdn.gene(1:200))')
            end
        end
    end
    
    answer = questdlg('Run enrichment analysis with top 200 up-regulated DE genes?');
    if strcmp(answer,'Yes')
        gui.i_enrichtest(Tup.gene(1:min(numel(Tup.gene),200)));
    end
    %pause(3);
    answer = questdlg('Run enrichment analysis with top 200 down-regulated DE genes?');
    if strcmp(answer,'Yes')
        gui.i_enrichtest(Tdn.gene(1:min(numel(Tdn.gene),200)));
    end

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