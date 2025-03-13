function callback_DEGene2Groups(src, ~)

    isatac = false;
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    % if ~gui.gui_showrefinfo('DE Analysis', FigureHandle), return; end

    %[wkdir] = gui.i_getwrkdir;
    %if isempty(wkdir), return; end
    extprogname = 'scgeatool_DEAnalysis';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    
    [i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
    if isscalar(i1) || isscalar(i2), return; end
    
    % --------
    a=sprintf('%s vs. %s',cL1{1}, cL2{1});
    b=sprintf('%s vs. %s',cL2{1}, cL1{1});
    answer = gui.myQuestdlg(FigureHandle, 'Which vs. which?','',{a,b},a);
    switch answer
        case a
        case b
            i3=i1; i1=i2; i2=i3;
            cL3=cL1; cL1=cL2; cL2=cL3;
        otherwise
            return;
    end
    % ----------
    
    methodtag = "ranksum";
   
    try
        switch methodtag
            case 'ranksum'
                T = sc_deg(sce.X(:, i1), sce.X(:, i2), sce.g, 1, ...
                    true, FigureHandle);
            case 'deseq2'
                [ok] = gui.i_confirmscript('DE analysis (DESeq2)', ...
                    'R_DESeq2', 'r');
                if ~ok, return; end
                fw = gui.myWaitbar(FigureHandle);
                T = run.r_DESeq2(sce.X(:, i1), sce.X(:, i2), sce.g);
                gui.myWaitbar(FigureHandle, fw);
            case 'mast'
                [ok] = gui.i_confirmscript('DE analysis (MAST)', 'R_MAST', 'r');
                if ~ok, return; end
                fw = gui.myWaitbar(FigureHandle);
                T = run.r_MAST(sce.X(:, i1), sce.X(:, i2), sce.g);
                gui.myWaitbar(FigureHandle, fw);
        end

    catch ME
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
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
        T = sortrows(T, 'p_val_adj', 'ascend');
        T = sortrows(T, 'pct_1', 'ascend');
        T = sortrows(T, 'pct_2', 'descend');
        T = sortrows(T, 'avg_log2FC', 'ascend');
    catch ME
        warning(ME.message);
    end

    try
        % avg_1 = mean(X,2);
        % avg_2 = mean(Y,2);
        % pct_1 = sum(X>0,2)./size(X,2);
        % pct_2 = sum(Y>0,2)./size(Y,2);
        if contains(T.Properties.VariableNames{5}, 'avg_1')
            T.Properties.VariableNames{5} = sprintf('%s_%s', ...
                T.Properties.VariableNames{5}, ...
                matlab.lang.makeValidName(string(cL1)));
        end

        if contains(T.Properties.VariableNames{6}, 'avg_2')
            T.Properties.VariableNames{6} = sprintf('%s_%s', ...
                T.Properties.VariableNames{6}, ...
                matlab.lang.makeValidName(string(cL2)));
        end

        if contains(T.Properties.VariableNames{7}, 'pct_1')
            T.Properties.VariableNames{7} = sprintf('%s_%s', ...
                T.Properties.VariableNames{7}, ...
                matlab.lang.makeValidName(string(cL1)));
        end

        if contains(T.Properties.VariableNames{8}, 'pct_2')
            T.Properties.VariableNames{8} = sprintf('%s_%s', ...
                T.Properties.VariableNames{8}, ...
                matlab.lang.makeValidName(string(cL2)));
        end
    catch ME
        warning(ME.message);
    end

    outfile = sprintf('%s_vs_%s_DE_results', ...
        matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
    if isatac, T.gene = "chr" + T.gene; end

    [filetype, filesaved] = gui.i_exporttable(T, true, ...
        'Tdegenelist', outfile, [], "All_genes", FigureHandle);

    tf = 0;
    if ~(ismcc || isdeployed) && strcmp(filetype, 'Workspace')
        [Tup, Tdn] = pkg.e_processDETable(T, [], FigureHandle);
        if isempty(Tup) && isempty(Tdn), return; end
        labels = {'Save DE results (selected up-regulated) to variable named:', ...
            'Save DE results (selected down-regulated) to variable named:'};
        vars = {'Tup', 'Tdn'};
        values = {Tup, Tdn};
        [~, tf] = export2wsdlg(labels, vars, values);
        if tf ~= 1, return; end
    end

    if ~isempty(filesaved)
        if strcmp(filetype, 'Excel file')
            %answer = gui.myQuestdlg(FigureHandle, 'Save up- and down-regulated genes to seperate sheets?');
            %if strcmp(answer, 'Yes')
                [Tup, Tdn] = pkg.e_processDETable(T,[],FigureHandle);
                % strcmp(extractAfter(filesaved,strlength(filesaved)-4),'xlsx')
                if ~isempty(Tup)
                    writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
                end
                if ~isempty(Tdn)
                    writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
                end
                gui.myHelpdlg(FigureHandle, sprintf('Result has been saved in %s', filesaved));
                %writetable(Tup,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet',);
                %writetable(Tdn,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet');
                tf = 1;
            %end
        elseif strcmp(filetype, 'Text file')
            % strcmp(extractAfter(filesaved,strlength(filesaved)-3),'txt')
            answer = gui.myQuestdlg(FigureHandle, 'Save up- and down-regulated genes to seperate files?');
            if strcmp(answer, 'Yes')
                [Tup, Tdn] = pkg.e_processDETable(T,[],FigureHandle);
                if ~isempty(Tup)
                    [~, ~] = gui.i_exporttable(Tup, true, 'Tup', 'Upregulated', 'Text file','', FigureHandle);
                end
                if ~isempty(Tdn)
                    [~, ~] = gui.i_exporttable(Tdn, true, 'Tdn', 'Downregulated', 'Text file','', FigureHandle);
                end
                tf = 1;
            end
        end
    end

    if tf ~= 1, return; end
    disp('To run enrichment analysis, type:');
    disp('run.web_Enrichr(Tup.gene(1:250))');
    disp('run.web_Enrichr(Tdn.gene(1:250))');

    [outgenelist, outbackgroundlist, enrichrtype] = ...
        gui.gui_prepenrichr(Tup.gene, sce.g,... 
       'Run enrichment analysis with up-regulated DE genes?', FigureHandle);

    if ~isempty(outbackgroundlist)
        gui.callback_RunEnrichr(src, [], outgenelist, enrichrtype, ...
            outbackgroundlist, "Up");
    end
    
    [outgenelist, outbackgroundlist, enrichrtype] = ...
        gui.gui_prepenrichr(Tdn.gene, sce.g,... 
       'Run enrichment analysis with down-regulated DE genes?', FigureHandle);

    if ~isempty(outbackgroundlist)
        gui.callback_RunEnrichr(src, [], outgenelist, enrichrtype, ...
            outbackgroundlist, "Down");
    end
    
end
