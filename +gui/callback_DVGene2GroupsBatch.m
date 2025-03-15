function callback_DVGene2GroupsBatch(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('DV in Batch Mode', FigureHandle), return; end
    
    extprogname = 'scgeatool_DVAnalysis_Batch';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end

prefixtag = 'DV';

a=sce.NumGenes;
[sce] = gui.i_selectinfogenes(sce, [], FigureHandle);
b=sce.NumGenes;
fprintf('%d genes removed.\n', a-b);

[done, CellTypeList, i1, i2, cL1, cL2, ...
    outdir] = gui.i_batchmodeprep(sce, prefixtag, wrkdir, FigureHandle);
if ~done, return; end

%[runenrichr] = gui.i_enrichrprep;
[runenrichr] = gui.myQuestdlg(FigureHandle, 'Run Enrichr with top 250 DV genes? Results will be saved in the output Excel files.','');
if strcmp(runenrichr,'Cancel'), return; end

fw = gui.myWaitbar(FigureHandle);
for k=1:length(CellTypeList)
   
    gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('Processing %s ...', CellTypeList{k}), ...
        (k-1)/length(CellTypeList));

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

    outfile = sprintf('%s_%s_vs_%s_%s.xlsx', ...
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
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
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
gui.myWaitbar(FigureHandle, fw);

answer=gui.myQuestdlg(FigureHandle, sprintf('Result files saved. Open the folder %s?', outdir), '');
if strcmp(answer,'Yes'), winopen(outdir); end

end