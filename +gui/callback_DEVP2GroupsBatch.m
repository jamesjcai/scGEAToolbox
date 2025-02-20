function callback_DEVP2GroupsBatch(src, ~)

[~, sce] = gui.gui_getfigsce(src);
%if ~gui.gui_showrefinfo('DE in Batch Mode'), return; end
extprogname = 'scgeatool_DEVPAnalysis_Batch';
preftagname = 'externalwrkpath';
[wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wrkdir), return; end

prefixtag = 'DEVP';
[done, CellTypeList, i1, i2, cL1, cL2, ... 
    outdir] = gui.i_batchmodeprep(sce, prefixtag, wrkdir);
if ~done, return; end

[paramset] = gui.i_degparamset(true);

% ------------------------------------------ DE
fw = gui.gui_waitbar_adv;
for k=1:length(CellTypeList)
   
    gui.gui_waitbar_adv(fw, ...
        (k-1)/length(CellTypeList), ...
        sprintf('Processing %s ...', CellTypeList{k}));

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
                   sce.g, 1, false);
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
        if length(Item) == length(Description)
            Tnt = table(Item, Description);
        else
            Tnt = table(Item);
        end
        
        [Tup, Tdn] = pkg.e_processDETable(T, paramset);
        try
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_genes');
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end
        
        try
            [Tlist1] = run.ml_Enrichr(Tup.gene(1:min([250 height(Tup)])), ...
                        T.gene, ["GO_Biological_Process_2023", ...
                                 "GO_Molecular_Function_2023"]);
            Tbp1 = Tlist1{1};
            Tmf1 = Tlist1{2};
            in_writetable(Tbp1, filesaved, 'Up_250_GO_BP');
            in_writetable(Tmf1, filesaved, 'Up_250_GO_MF');
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
end
gui.gui_waitbar_adv(fw);

% ------------------------------------------ DV
fw = gui.gui_waitbar_adv;
for k=1:length(CellTypeList)
    gui.gui_waitbar_adv(fw, ...
        (k-1)/length(CellTypeList), ...
        sprintf('Processing %s ...', CellTypeList{k}));
    idx = sce.c_cell_type_tx == CellTypeList{k};
    sce1 = sce.selectcells(i1&idx);
    sce1 = sce1.qcfilter;

    sce2 = sce.selectcells(i2&idx);
    sce2 = sce2.qcfilter;

    if sce1.NumCells < 10 || sce2.NumCells < 10 || sce1.NumGenes < 10 || sce2.NumGenes < 10
        warning('Filtered SCE contains too few cells (n < 10) or genes (n < 10).');
        continue;
    end

    [T] = gui.e_dvanalysis_splinefit(sce1, sce2, cL1, cL2);
    outfile = sprintf('%s_DV_%s_vs_%s_%s.xlsx', ...
        prefixtag,...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)), ...
        matlab.lang.makeValidName(string(CellTypeList{k})));
        filesaved = fullfile(outdir, outfile);        
         
        Tup = T(T.DiffSign > 0, :);
        Tdn = T(T.DiffSign < 0, :);

        Item = T.Properties.VariableNames';
        Item = [Item; {'# of cells in sample 1';'# of cells in sample 2'}];
        
        Description = {'gene name';'log mean in sample 1';...
            'log CV in sample 1'; 'dropout rate in sample 1';...
            'distance to curve 1';'p-value of distance in sample 1';...
            'FDR of distance in sample 1';'log mean in sample 2';...
            'log CV in sample 2'; 'dropout rate in sample 2';...
            'distance to curve 2'; 'p-value of distance in sample 2';...
            'FDR of distance in sample 2'; 'Difference in distances';...
            'Sign of difference';'p-value of DV test';...
            sprintf('%d',sce1.NumCells); sprintf('%d',sce2.NumCells)};
        if length(Item) == length(Description)
            Tnt = table(Item, Description);
        else
            assignin("base","Item", Item);
            assignin("base","Description", Description);
            Tnt = table(Item);
            warning('Variables must have the same number of rows.');
        end

        try
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end

        try

            [Tlist1] = run.ml_Enrichr(Tup.gene(1:min([250 height(Tup)])), ...
                        T.gene, ["GO_Biological_Process_2023", ...
                                 "GO_Molecular_Function_2023",...
                                 "Reactome_Pathways_2024",...
                                 "KEGG_2021_Human"]);
                                 % KEGG_2019_Mouse
            Tbp1 = Tlist1{1};
            Tmf1 = Tlist1{2};
            in_writetable(Tbp1, filesaved, 'Up_250_GO_BP');
            in_writetable(Tmf1, filesaved, 'Up_250_GO_MF');
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
gui.gui_waitbar_adv(fw);

% ----------------------------- DP
ctag = {"H", "C2", "C5", "C6", "C7"};
pw1 = fileparts(mfilename('fullpath'));

for c = 1:length(ctag)
    dbfile = fullfile(pw1, '..', 'resources', 'MSigDB', ...
                    sprintf('msigdb_%s.mat', ctag{c}));
    load(dbfile,'setmatrx','setnames','setgenes');
    
    fw = gui.gui_waitbar_adv;
    sceX = log1p(sc_norm(sce.X));
    for k=1:length(CellTypeList)
        gui.gui_waitbar_adv(fw, ...
            (k-1)/length(CellTypeList), ...
            sprintf('Processing %s ...', CellTypeList{k}));
    
        outfile = sprintf('%s_DP_%s_%s_vs_%s_%s.xlsx', ...
            prefixtag, ctag{c}, ...
            matlab.lang.makeValidName(string(cL1)), ...
            matlab.lang.makeValidName(string(cL2)), ...
            matlab.lang.makeValidName(string(CellTypeList{k})));
            filesaved = fullfile(outdir, outfile);
    
        idx = sce.c_cell_type_tx == CellTypeList{k};
        try
            T = sc_dpg(sceX(:, i1&idx), sceX(:, i2&idx), sce.g, ...
                setmatrx, setnames, setgenes);
            Tup = T(T.avg_log2FC > 0, :);
            Tdn = T(T.avg_log2FC < 0, :);
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All programs');
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
        catch ME
            warning(ME.message);
        end
    end
    gui.gui_waitbar_adv(fw);
end

answer=questdlg(sprintf('Result files saved. Open the folder %s?', outdir), '');
if strcmp(answer,'Yes'), winopen(outdir); end
end   % end of function


function in_writetable(Tmf1, filesaved, shtname)
    if ~isempty(Tmf1) && istable(Tmf1) && height(Tmf1) > 0
        writetable(Tmf1, filesaved, "FileType", "spreadsheet", 'Sheet', shtname);
    end
end

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