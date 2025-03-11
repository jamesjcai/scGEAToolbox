function callback_DEGene2GroupsBatch(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('DE in Batch Mode'), return; end

    extprogname = 'scgeatool_DEAnalysis_Batch';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
    if isempty(wrkdir), return; end

prefixtag = 'DE';
[done, CellTypeList, i1, i2, cL1, cL2,... 
    outdir] = gui.i_batchmodeprep(sce, prefixtag, wrkdir, FigureHandle);
if ~done, return; end

%[runenrichr] = gui.i_enrichrprep;
[runenrichr] = gui.myQuestdlg(FigureHandle, 'Run Enrichr with top 250 DE genes? Results will be saved in the output Excel files.','');
if strcmp(runenrichr,'Cancel'), return; end

[paramset] = gui.i_degparamset;

fw = gui.gui_waitbar_adv;
for k=1:length(CellTypeList)
   
    gui.gui_waitbar_adv(fw, ...
        (k-1)/length(CellTypeList), ...
        sprintf('Processing %s ...', CellTypeList{k}));

    outfile = sprintf('%s_%s_vs_%s_%s.xlsx', ...
        prefixtag, ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)), ...
        matlab.lang.makeValidName(string(CellTypeList{k})));
        filesaved = fullfile(outdir, outfile);

    idx = sce.c_cell_type_tx == CellTypeList{k};
    T = [];
    try
        T = sc_deg(sce.X(:, i1&idx), sce.X(:, i2&idx), ...
                   sce.g, 1, false, FigureHandle);
    catch ME
        disp(ME.message);
    end

    if ~isempty(T)
        T = in_DETableProcess(T, cL1, cL2);

        Item = T.Properties.VariableNames';
        Item = [Item; {'# of cells in sample 1';'# of cells in sample 2'}];
        Description = {'gene name';'p-value';...
            'log2 fold change between average expression';...
            'absolute value of log2 fold change';...
            'average expression in sample 1';...
            'average expression in sample 2';...
            'percentage of cells expressing the gene in sample 1';...
            'percentage of cells expressing the gene in sample 2';...
            'adjusted p-value'; 'test statistic'; ...
             sprintf('%d',sum(i1&idx)); sprintf('%d',sum(i2&idx))};
        %assignin('base','Item', Item);
        %assignin('base','Description', Description);        
        Tnt = table(Item, Description);
        
        [Tup, Tdn] = pkg.e_processDETable(T, paramset, FigureHandle);
        try
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_genes');
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end
        
        % - start of enrichr
        if ~isempty(runenrichr) && strcmp(runenrichr, 'Yes')
            try
                % [Tbp1, Tmf1]= run.r_enrichR(Tup.gene(1:min([250 height(Tup)])));
                %[Tbp1, Tmf1] = run.py_GSEApy_enr(Tup.gene(1:min([250 height(Tup)])), ...
                %    T.gene, tempdir);

                [Tlist1] = run.ml_Enrichr(Tup.gene(1:min([250 height(Tup)])), ...
                            T.gene, ["GO_Biological_Process_2023", ...
                                     "GO_Molecular_Function_2023"]);
                Tbp1 = Tlist1{1};
                Tmf1 = Tlist1{2};
                in_writetable(Tbp1, filesaved, 'Up_250_GO_BP');
                in_writetable(Tmf1, filesaved, 'Up_250_GO_MF');
                % [Tbp2, Tmf2] = run.r_enrichR(Tdn.gene(1:min([250 height(Tdn)])));
                % [Tbp2, Tmf2] = run.py_GSEApy_enr(Tdn.gene(1:min([250 height(Tdn)])), ...
                %    T.gene, tempdir);

                [Tlist2] = run.ml_Enrichr(Tdn.gene(1:min([250 height(Tdn)])), ...
                            T.gene, ["GO_Biological_Process_2023", ...
                                     "GO_Molecular_Function_2023"]);
                Tbp2 = Tlist2{1};
                Tmf2 = Tlist2{2};
                in_writetable(Tbp2, filesaved, 'Dn_250_GO_BP');
                in_writetable(Tmf2, filesaved, 'Dn_250_GO_MF');
            catch ME
                warning(ME.message);
            end
        end
        % - end of enrichr

    end
end
gui.gui_waitbar_adv(fw);

answer=gui.myQuestdlg(FigureHandle, sprintf('Result files saved. Open the folder %s?', outdir), '');
if strcmp(answer,'Yes'), winopen(outdir); end

    function in_writetable(Tmf1, filesaved, shtname)
        if ~isempty(Tmf1) && istable(Tmf1) && height(Tmf1) > 0
            writetable(Tmf1, filesaved, "FileType", "spreadsheet", 'Sheet', shtname);
        end
    end
end   % end of function



function [T] = in_DETableProcess(T, cL1, cL2)
    try
        T = sortrows(T, 'p_val_adj', 'ascend');
        T = sortrows(T, 'pct_1', 'ascend');
        T = sortrows(T, 'pct_2', 'descend');
        T = sortrows(T, 'avg_log2FC', 'ascend');
    catch ME
        warning(ME.message);
    end
    
    try
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
        variables = T.Properties.VariableNames;
        for k = 1:length(variables)
            xx = T.(variables{k});
            if isnumeric(xx) && any(isinf(xx))
                xx(isinf(xx) & xx > 0) = 1e99;
                xx(isinf(xx) & xx < 0) = -1e99;
                T.(variables{k}) = xx;
            end
        end
end