function callback_DEVP2GroupsBatch(src, ~)

    [FigureHandle, sce_ori] = gui.gui_getfigsce(src);
    sce = copy(sce_ori);

    extprogname = 'scgeatool_DEVPAnalysis_Batch';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end
    
    prefixtag = 'DEVP';
    [done, CellTypeList, i1, i2, cL1, cL2, ... 
        outdir] = gui.i_batchmodeprep(sce, prefixtag, ...
                wrkdir, FigureHandle);
    if ~done, return; end
    
    answer = gui.myQuestdlg(FigureHandle, "Set DE gene filter parameters?", ...
        "DE Genes", {'Yes','No, use previous','Cancel'}, 'Yes');
    if isempty(answer) || strcmp(answer, 'Cancel'), return; end
    if strcmp(answer, 'Yes')
        [paramset] = gui.i_degparamset(false, FigureHandle);
    else
        [paramset] = gui.i_degparamset(true, FigureHandle);
    end
    
    % ------------------------------------------ DE
    fw = gui.myWaitbar(FigureHandle);
    for k=1:length(CellTypeList)
       
        gui.myWaitbar(FigureHandle, fw, false, '', ...
            sprintf('DE - Processing %s ...', CellTypeList{k}), ...
            (k-0.5)/length(CellTypeList));
    
        outfile = sprintf('%s_DE_%s_vs_%s_%s.xlsx', ...
            prefixtag, ...
            matlab.lang.makeValidName(string(cL1)), ...
            matlab.lang.makeValidName(string(cL2)), ...
            matlab.lang.makeValidName(string(CellTypeList{k})));
            filesaved = fullfile(outdir, outfile);
    
        idx = sce.c_cell_type_tx == CellTypeList{k};
        T = [];
        try
            T = sc_deg(sce.X(:, i1&idx), ...
                       sce.X(:, i2&idx), ...
                       sce.g, 1, false, FigureHandle);
        catch ME
            disp(ME.message);
        end
    
        if ~isempty(T)
            [T, Tnt] = pkg.in_DETableProcess(T, cL1, cL2, sum(i1&idx), sum(i2&idx));
            
            % Item = T.Properties.VariableNames';
            % Item = [Item; {'# of cells in sample 1';'# of cells in sample 2'}];
            % Description = {'gene name';'p-value';...
            %     'log2 fold change between average expression';...
            %     'absolute value of log2 fold change';...
            %     'average expression in sample 1';...
            %     'average expression in sample 2';...
            %     'percentage of cells expressing the gene in sample 1';...
            %     'percentage of cells expressing the gene in sample 2';...
            %     'adjusted p-value'; 'test statistic'; ...
            %      sprintf('%d',sum(i1&idx)); sprintf('%d',sum(i2&idx))};
            % if length(Item) == length(Description)
            %     Tnt = table(Item, Description);
            % else
            %     Tnt = table(Item);
            % end
            
            [Tup, Tdn] = pkg.e_processdetable(T, paramset, FigureHandle);
            try
                gui.e_tupdn2xlsx(Tup, Tdn, T, filesaved);
                writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
                gui.e_enrichrxlsx(Tup, Tdn, T, filesaved);
            catch ME
                warning(ME.message);
            end
        end
        
        % try
        %     [done] = llm.e_DETableSummarizer(Tbp1,Tmf1,Tbp2,Tmf2, ...
        %         sprintf('DE_%s', CellTypeList{k}));
        % catch
        % 
        % end
        
    end
%   gui.myWaitbar(FigureHandle, fw);    
% ------------------------------------------ DV
%   fw = gui.myWaitbar(FigureHandle);
    for k=1:length(CellTypeList)
        gui.myWaitbar(FigureHandle, fw, false, '', ...
            sprintf('DV - Processing %s ...', CellTypeList{k}), ...
            (k-0.5)/length(CellTypeList));
        idx = sce.c_cell_type_tx == CellTypeList{k};
        sce1=copy(sce);
        sce1 = sce1.selectcells(i1&idx); % OK
        sce1 = sce1.qcfilter; % OK
    
        sce2 = copy(sce);
        sce2 = sce2.selectcells(i2&idx); % OK
        sce2 = sce2.qcfilter; % OK
        
        notok = false;
        if sce1.NumCells < 10
            disp(CellTypeList{k})
            warning('Filtered SCE 1 contains too few cells (NumCells < 10)');
            notok = true;
        end
        if  sce2.NumCells < 10
            disp(CellTypeList{k})
            warning('Filtered SCE 2 contains too few cells (NumCells < 10)');
            notok = true;
        end
        if sce1.NumGenes < 10
            disp(CellTypeList{k})
            warning('Filtered SCE 1 contains too few genes (NumGenes < 10)');
            notok = true;
        end
        if sce2.NumGenes < 10
            disp(CellTypeList{k})
            warning('Filtered SCE 2 contains too few genes (NumGenes < 10)');
            notok = true;
        end
        if notok, continue; end
    
        [T] = gui.e_dvanalysis_splinefit(sce1, sce2, cL1, cL2);

        outfile = sprintf('%s_DV_%s_vs_%s_%s.xlsx', ...
            prefixtag,...
            matlab.lang.makeValidName(string(cL1)), ...
            matlab.lang.makeValidName(string(cL2)), ...
            matlab.lang.makeValidName(string(CellTypeList{k})));
            filesaved = fullfile(outdir, outfile);        
             
            Tup = T(T.DiffSign > 0, :);
            Tdn = T(T.DiffSign < 0, :);
            
            [T, Tnt] = pkg.in_DVTableProcess(T, cL1, cL2);

            % Item = T.Properties.VariableNames';
            % Item = [Item; {'# of cells in sample 1';'# of cells in sample 2'}];
            % 
            % Description = {'gene name';'log mean in sample 1';...
            %     'log CV in sample 1'; 'dropout rate in sample 1';...
            %     'distance to curve 1';'p-value of distance in sample 1';...
            %     'FDR of distance in sample 1';'log mean in sample 2';...
            %     'log CV in sample 2'; 'dropout rate in sample 2';...
            %     'distance to curve 2'; 'p-value of distance in sample 2';...
            %     'FDR of distance in sample 2'; 'Difference in distances';...
            %     'Sign of difference';'p-value of DV test';...
            %     sprintf('%d',sce1.NumCells); sprintf('%d',sce2.NumCells)};
            % if length(Item) == length(Description)
            %     Tnt = table(Item, Description);
            % else
            %     assignin("base","Item", Item);
            %     assignin("base","Description", Description);
            %     Tnt = table(Item);
            %     warning('Variables must have the same number of rows.');
            % end
    
            try
                gui.e_tupdn2xlsx(Tup,Tdn,T,filesaved);
                writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');                
                gui.e_enrichrxlsx(Tup,Tdn,T,filesaved);
            catch ME
                warning(ME.message);
            end
    end
    % gui.myWaitbar(FigureHandle, fw);
    
    % ----------------------------- DP
    ctag = {"H", "C2", "C5", "C6", "C7"}';
    ccat = {"H: Hallmark gene sets (broadly defined, high-quality gene signatures representing specific biological states or processes)", ...
        "C2: Curated gene sets (pathways from KEGG, Reactome, BioCarta, and literature)",...
        "C5: Gene Ontology (GO) gene sets (BP: biological process, CC: cellular component, MF: molecular function)",...
        "C6: Oncogenic signatures (gene sets linked to cancer-related mutations and pathways)",...
        "C7: Immunologic signatures (gene sets related to immune cell expression and responses)"}';

    pw1 = fileparts(mfilename('fullpath'));
    
    for c = 1:length(ctag)
        dbfile = fullfile(pw1, '..', 'assets', 'MSigDB', ...
                        sprintf('msigdb_%s.mat', ctag{c}));
        load(dbfile,'setmatrx','setnames','setgenes');
        
%       fw = gui.myWaitbar(FigureHandle);
        sceX = log1p(sc_norm(sce.X));
        for k=1:length(CellTypeList)
            gui.myWaitbar(FigureHandle, fw, false, '', ...
                sprintf('DP - Processing %s ...', CellTypeList{k}), ...
                (k-0.5)/length(CellTypeList));
        
            outfile = sprintf('%s_DP_%s_vs_%s_%s.xlsx', ...
                prefixtag, ...
                matlab.lang.makeValidName(string(cL1)), ...
                matlab.lang.makeValidName(string(cL2)), ...
                matlab.lang.makeValidName(string(CellTypeList{k})));
                filesaved = fullfile(outdir, outfile);
        
            idx = sce.c_cell_type_tx == CellTypeList{k};
            try
                T = sc_dpg(sceX(:, i1&idx), sceX(:, i2&idx), sce.g, ...
                    setmatrx, setnames, setgenes);
                if ~isempty(T)
                    writetable(T, filesaved, 'FileType', 'spreadsheet', ...
                        'Sheet', sprintf('All_%s', ctag{c}));
                    Tup = T(T.avg_log2FC > 0, :);                    
                    Tdn = T(T.avg_log2FC < 0, :);
                    if ~isempty(Tup)
                        writetable(Tup, filesaved, "FileType", "spreadsheet", ...
                            'Sheet', sprintf('Up_%s', ctag{c}));
                    end
                    if ~isempty(Tdn)
                        writetable(Tdn, filesaved, "FileType", "spreadsheet", ...
                            'Sheet', sprintf('Dn_%s', ctag{c}));
                    end
                end
            catch ME
                warning(ME.message);
            end
        end
        % gui.myWaitbar(FigureHandle, fw);
    end    

    Tnt = table(ctag, ccat);
    writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');

    gui.myWaitbar(FigureHandle, fw);
    
    % answer = gui.myQuestdlg(FigureHandle, 'Use LLM to generate enrichment analysis report?', '');
    % if strcmp(answer,'Yes')
    %     gui.sc_llm_enrichr2word(outdir);
    % end
    % 
    % answer = gui.myQuestdlg(FigureHandle, sprintf('Result files saved. Open the folder %s?', outdir), '');
    % if strcmp(answer,'Yes'), winopen(outdir); end

    items = {'LLM Summarize', 'Open Output Folder'};
    selected = gui.myChecklistdlg(FigureHandle, items, ...
        'Title', 'Select Items','DefaultSelection', [1 2]);
    if isempty(selected), return; end
    
    if any(contains(selected, 'LLM Summarize'))
        %fw = gui.myWaitbar(FigureHandle, [], false, 'Use LLM to generate enrichment analysis report');        
        gui.sc_llm_enrichr2word(outdir);
        %gui.myWaitbar(FigureHandle, fw);
    end

    if any(contains(selected, 'Open Output Folder'))        
        if strcmp('Yes', gui.myQuestdlg(FigureHandle,'Open Output Folder?'))            
            winopen(outdir);
        end
    end    

end

