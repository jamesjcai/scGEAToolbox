function callback_DPGene2GroupsBatch(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('DP in Batch Mode', FigureHandle), return; end

extprogname = 'scgeatool_DPAnalysis_Batch';
preftagname = 'externalwrkpath';
[wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wrkdir), return; end

prefixtag = 'DP';

[sce, spciestag] = gui.i_selectinfogenes(sce, [], FigureHandle);
if isempty(spciestag), return; end
[done, CellTypeList, i1, i2, cL1, cL2,...
outdir] = gui.i_batchmodeprep(sce, prefixtag, wrkdir, FigureHandle);
if ~done, return; end

[indx1, species] = gui.i_selgenecollection(FigureHandle, spciestag);
if isempty(indx1), return; end
[setmatrx, setnames, setgenes] = pkg.e_getgenesets(indx1, species, FigureHandle); % (indx1);
if isempty(setmatrx) || isempty(setnames) || isempty(setgenes)
    return;
end

fw = gui.myWaitbar(FigureHandle);
minCellsPerGroup = 2;
allSaved = true;
anySaved = false;
numSkipped = 0;
numFailed = 0;

ranknorm   = true;
bgsubtract = true;
sceX = log1p(sc_norm(sce.X));
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
    n1 = sum(i1 & idx);
    n2 = sum(i2 & idx);
    if n1 < minCellsPerGroup || n2 < minCellsPerGroup
        warning(['Skipping DP batch analysis for %s: sample sizes are too ' ...
            'small (%s=%d, %s=%d; minimum=%d).'], CellTypeList{k}, ...
            string(cL1), n1, string(cL2), n2, minCellsPerGroup);
        allSaved = false;
        numSkipped = numSkipped + 1;
        continue;
    end
    try
        T = sc_dpg(sceX(:, i1&idx), sceX(:, i2&idx), sce.g, ...
            setmatrx, setnames, setgenes, ranknorm, bgsubtract);
        Tup = T(T.avg_log2FC > 0, :);
        Tdn = T(T.avg_log2FC < 0, :);
        writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All programs');
        writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
        writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
        anySaved = true;
    catch ME
        warning(ME.message);
        allSaved = false;
        numFailed = numFailed + 1;
    end
end
gui.myWaitbar(FigureHandle, fw);

if ~anySaved
    gui.myHelpdlg(FigureHandle, ...
        sprintf(['No DP batch result files were saved.\n' ...
        'Skipped cell types: %d\nFailed cell types: %d'], ...
        numSkipped, numFailed));
    return;
end

if allSaved
    msg = sprintf('Result files saved. Open the folder %s?', outdir);
else
    msg = sprintf(['Some DP batch result files were not saved.\n' ...
        'Skipped cell types: %d\nFailed cell types: %d\n' ...
        'Open the folder %s?'], numSkipped, numFailed, outdir);
end
answer=gui.myQuestdlg(FigureHandle, msg, '');
if strcmp(answer,'Yes'), winopen(outdir); end

answer = gui.myQuestdlg(FigureHandle, 'Use LLM to generate program analysis report?', '');
if strcmp(answer,'Yes')
    gui.sc_llm_dp2word(outdir, FigureHandle);
end

end
