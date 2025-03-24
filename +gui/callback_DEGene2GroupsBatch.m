function callback_DEGene2GroupsBatch(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('DE in Batch Mode', FigureHandle), return; end

    extprogname = 'scgeatool_DEAnalysis_Batch';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end

prefixtag = 'DE';
[done, CellTypeList, i1, i2, cL1, cL2,... 
    outdir] = gui.i_batchmodeprep(sce, prefixtag, wrkdir, FigureHandle);
if ~done, return; end

%[runenrichr] = gui.i_enrichrprep;
[runenrichr] = gui.myQuestdlg(FigureHandle, ...
    ['Run Enrichr with top 250 DE genes? Results will ' ...
    'be saved in the output Excel files.'],'');
if strcmp(runenrichr,'Cancel'), return; end

[paramset] = gui.i_degparamset(false, FigureHandle);

fw = gui.myWaitbar(FigureHandle);
for k=1:length(CellTypeList)
   
    gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('Processing %s ...', CellTypeList{k}), ...
        (k-1)/length(CellTypeList));

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
        [T, Tnt] = pkg.in_DETableProcess(T, cL1, cL2, sum(i1&idx), sum(i2&idx));
        
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
                gui.e_enrichrxlsx(Tup,Tdn,T,filesaved);

            catch ME
                warning(ME.message);
            end
        end
        % - end of enrichr
    end
end
gui.myWaitbar(FigureHandle, fw);

answer=gui.myQuestdlg(FigureHandle, sprintf('Result files saved. Open the folder %s?', outdir), '');
if strcmp(answer,'Yes'), winopen(outdir); end

    % function in_writetable(Tmf1, filesaved, shtname)
    %     if ~isempty(Tmf1) && istable(Tmf1) && height(Tmf1) > 0
    %         writetable(Tmf1, filesaved, "FileType", "spreadsheet", 'Sheet', shtname);
    %     end
    % end

end   % end of function





