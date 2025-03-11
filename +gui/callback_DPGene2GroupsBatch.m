function callback_DPGene2GroupsBatch(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    if ~gui.gui_showrefinfo('DP in Batch Mode', FigureHandle), return; end

    extprogname = 'scgeatool_DPAnalysis_Batch';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end

    prefixtag = 'DP';
    
    [sce] = gui.i_selectinfogenes(sce, [], FigureHandle);
    [done, CellTypeList, i1, i2, cL1, cL2,... 
        outdir] = gui.i_batchmodeprep(sce, prefixtag, wrkdir);
    if ~done, return; end
    
    [indx1, species] = gui.i_selgenecollection(FigureHandle);
    if isempty(indx1), return; end
    [setmatrx, setnames, setgenes] = pkg.e_getgenesets(indx1, species); %(indx1);
    if isempty(setmatrx) || isempty(setnames) || isempty(setgenes) 
        return; 
    end
    
    fw = gui.gui_waitbar_adv;
    
    sceX = log1p(sc_norm(sce.X));
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
    
    answer=gui.myQuestdlg(FigureHandle, sprintf('Result files saved. Open the folder %s?', outdir), '');
    if strcmp(answer,'Yes'), winopen(outdir); end
end