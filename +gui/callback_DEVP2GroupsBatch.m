function callback_DEVP2GroupsBatch(src, ~)

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

extprogname = 'scgeatool_DEVPAnalysis_Batch';
preftagname = 'externalwrkpath';
[wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wrkdir), return; end

prefixtag = 'DEVP';
% The three infixes this callback writes. Without them the overwrite
% check in i_batchmodeprep probes DEVP_<g1>_vs_<g2>_<ct>.xlsx, which this
% callback never writes, so it could not fire.
[done, CellTypeList, i1, i2, cL1, cL2, ...
outdir] = gui.i_batchmodeprep(sce, prefixtag, ...
            wrkdir, FigureHandle, {'_DE', '_DV', '_DP'});
if ~done, return; end

% Enrichr is a web service, so it is asked for, as the DE and DV batch
% callbacks ask. It used to run in both phases unconditionally.
runenrichr = gui.myQuestdlg(FigureHandle, ...
    ['Run Enrichr with top 250 DE and DV genes? Results will ' ...
     'be saved in the output Excel files.'], '');
if ~ismember(runenrichr, {'Yes', 'No'}), return; end   % Cancel or closed
runenrichr = strcmp(runenrichr, 'Yes');

answer = gui.myQuestdlg(FigureHandle, "Set DE gene filter parameters?", ...
"DE Genes", {'Yes','No, use previous','Cancel'}, 'Yes');
if isempty(answer) || strcmp(answer, 'Cancel'), return; end
if strcmp(answer, 'Yes')
    [paramset] = gui.i_degparamset(false, FigureHandle);
else
    [paramset] = gui.i_degparamset(true, FigureHandle);
end
if isempty(paramset), return; end

% How the DV phase splits up from down. The file names stay DEVP_DV_...
% either way; the Note sheet records which was used.
direction = gui.i_dvdirection(FigureHandle);
if isempty(direction), return; end
% Likewise the DV p-value source: the Note sheet says which was used.
numPerm = gui.i_dvpermutations(FigureHandle, 'splinefit');
if isempty(numPerm), return; end

% MSigDB collections for the DP phase; needed here to size the progress bar.
ctag = {"H", "C2", "C5", "C6", "C7"}';

% One progress bar for the whole run. Each phase used to compute its own
% (k-0.5)/nCellTypes, so the bar started over for DE, DV and each of the
% five DP collections and ran backwards seven times in one run.
nCellType = length(CellTypeList);
nStepTotal = nCellType*(2 + numel(ctag));
iStep = 0;

% ------------------------------------------ DE
fw = gui.myWaitbar(FigureHandle);
% Closed on every exit: a failure inside the loop below left it open.
closeFw = onCleanup(@() gui.myWaitbar(FigureHandle, fw, true));
for k=1:length(CellTypeList)

    iStep = iStep + 1;
    gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('DE - Processing %s ...', CellTypeList{k}), ...
        (iStep - 0.5)/nStepTotal);

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

        [Tup, Tdn, ~, usedset] = pkg.e_processdetable(T, paramset, FigureHandle);
        Tnt = pkg.i_decutoffnote(Tnt, usedset);
        try
            gui.e_tupdn2xlsx(Tup, Tdn, T, filesaved);
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end
        if runenrichr
            try
                gui.e_enrichrxlsx(Tup, Tdn, T, filesaved);
            catch ME
                warning(ME.message);
            end
        end
    end

end
%   gui.myWaitbar(FigureHandle, fw);
% ------------------------------------------ DV
%   fw = gui.myWaitbar(FigureHandle);
for k=1:length(CellTypeList)
    iStep = iStep + 1;
    gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('DV - Processing %s ...', CellTypeList{k}), ...
        (iStep - 0.5)/nStepTotal);
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

    [T] = sc_dvg(sce1, sce2, cL1, cL2, 'splinefit', direction, NumPermutations=numPerm);

    outfile = sprintf('%s_DV_%s_vs_%s_%s.xlsx', ...
        prefixtag,...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)), ...
        matlab.lang.makeValidName(string(CellTypeList{k})));
        filesaved = fullfile(outdir, outfile);

        % A significance cutoff, as the DE branch above gets from
        % pkg.e_processdetable. Splitting on DiffSign alone put EVERY tested
        % gene into one list or the other - measured on this project's own
        % output, 7,550 of 7,550 - because sc_dvg returns the full ranked
        % table and nothing here ever filtered it. DiffDist > 0 additionally
        % drops the spline-boundary genes whose distance was discarded.
        % PKG.E_FDR has one output, so the trap that this comment used
        % to warn about -- E_FDR_BH returns h first and the adjusted p
        % fourth, and taking output 1 as a q-value silently inverts the
        % test -- cannot be sprung.
        dvq = pkg.e_fdr(T.pval);
        isok = T.DiffDist > 0 & dvq(:) <= 0.05;
        fprintf(['\nDV genes with BH q <= %.3f and a usable spline ' ...
            'distance are retained: %d of %d.\n'], 0.05, sum(isok), height(T));
        Tup = T(T.DiffSign > 0 & isok, :);
        Tdn = T(T.DiffSign < 0 & isok, :);

        [T, Tnt] = pkg.in_DVTableProcess(T, cL1, cL2, direction, numPerm);

        try
            gui.e_tupdn2xlsx(Tup,Tdn,T,filesaved);
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end
        if runenrichr
            try
                gui.e_enrichrxlsx(Tup,Tdn,T,filesaved);
            catch ME
                warning(ME.message);
            end
        end
end

% ----------------------------- DP
ccat = {"H: Hallmark gene sets (broadly defined, high-quality gene signatures representing specific biological states or processes)", ...
"C2: Curated gene sets (pathways from KEGG, Reactome, BioCarta, and literature)",...
"C5: Gene Ontology (GO) gene sets (BP: biological process, CC: cellular component, MF: molecular function)",...
"C6: Oncogenic signatures (gene sets linked to cancer-related mutations and pathways)",...
"C7: Immunologic signatures (gene sets related to immune cell expression and responses)"}';

pw1 = fileparts(mfilename('fullpath'));

ranknorm   = true;
bgsubtract = true;
sceX = log1p(sc_norm(sce.X));
for c = 1:length(ctag)
    dbfile = fullfile(pw1, '..', 'assets', 'MSigDB', ...
                    sprintf('msigdb_%s.mat', ctag{c}));
    load(dbfile,'setmatrx','setnames','setgenes');

    for k=1:length(CellTypeList)
        iStep = iStep + 1;
        gui.myWaitbar(FigureHandle, fw, false, '', ...
            sprintf('DP (%s) - Processing %s ...', ctag{c}, CellTypeList{k}), ...
            (iStep - 0.5)/nStepTotal);

        outfile = sprintf('%s_DP_%s_vs_%s_%s.xlsx', ...
            prefixtag, ...
            matlab.lang.makeValidName(string(cL1)), ...
            matlab.lang.makeValidName(string(cL2)), ...
            matlab.lang.makeValidName(string(CellTypeList{k})));
            filesaved = fullfile(outdir, outfile);

        idx = sce.c_cell_type_tx == CellTypeList{k};
        try
            T = sc_dpg(sceX(:, i1&idx), sceX(:, i2&idx), sce.g, ...
                setmatrx, setnames, setgenes, ranknorm, bgsubtract);
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
end

% The Note sheet goes in every DP file. It was written once, after the
% loops, into FILESAVED -- which by then named only the last cell type's
% file, so every other DP workbook went without it.
Tnt = table(ctag, ccat);
for k = 1:length(CellTypeList)
    outfile = sprintf('%s_DP_%s_vs_%s_%s.xlsx', ...
        prefixtag, ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)), ...
        matlab.lang.makeValidName(string(CellTypeList{k})));
    filesaved = fullfile(outdir, outfile);
    if isfile(filesaved)
        try
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end
    end
end

gui.myWaitbar(FigureHandle, fw);

% answer = gui.myQuestdlg(FigureHandle, 'Use LLM to generate enrichment analysis report?', '');
% if strcmp(answer,'Yes')
%     gui.sc_llm_enrichr2word(outdir);
% end
%
% answer = gui.myQuestdlg(FigureHandle, sprintf('Result files saved. Open the folder %s?', outdir), '');
% if strcmp(answer,'Yes'), winopen(outdir); end

% The LLM summary reads the Enrichr sheets, so it is offered only when
% Enrichr ran.
if runenrichr
    items = {'LLM Summarize', 'Open Output Folder'};
else
    items = {'Open Output Folder'};
end
selected = gui.myChecklistdlg(FigureHandle, items, ...
'Title', 'Select Items','DefaultSelection', 1:numel(items));
if isempty(selected), return; end

if any(contains(selected, 'LLM Summarize'))
    gui.sc_llm_enrichr2word(outdir, FigureHandle);
end

if any(contains(selected, 'Open Output Folder'))
    if strcmp('Yes', gui.myQuestdlg(FigureHandle,'Open Output Folder?'))
        winopen(outdir);
    end
end

end
