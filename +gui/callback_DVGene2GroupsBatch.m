function callback_DVGene2GroupsBatch(src, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
if ~gui.gui_showrefinfo('DV in Batch Mode'), return; end
prefixtag = 'DV';

a=sce.NumGenes;
[sce] = gui.i_selectinfogenes(sce);
b=sce.NumGenes;
fprintf('%d genes removed.\n', a-b);

[done, CellTypeList, i1, i2, cL1, cL2, ...
    outdir] = gui.i_batchmodeprep(sce, prefixtag);
if ~done, return; end

runenrichr = questdlg('Run Enrichr (R required) with top 250 DV genes? Results will be saved in the output Excel files.','');
if strcmp(runenrichr,'Cancel'), return; end

isok = false;
if strcmp(runenrichr, 'Yes')
    [isok] = gui.i_commoncheck_R('r_enrichR');
    if ~isok
        answer = questdlg('R Environment Error: It seems that your R environment is not configured correctly. This means that Gene Function Enrichment Analysis using enrichR cannot be performed for differentially expressed genes. Continue withouth enrichR?',''); 
        if ~strcmp(answer,'Yes'), return; end
    end
end

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
%{
    if sce1.NumCells < 50 || sce2.NumCells < 50
        warndlg('One of groups contains too few cells (n < 50). The result may not be reliable.','','modal');
    end
    if sce1.NumGenes < 50 || sce2.NumGenes < 50
        warndlg('One of groups contains too few genes (n < 50). The result may not be reliable.','','modal');
    end

    if ~isequal(sce1.g, sce2.g)
        [g_ori, ia, ib] = intersect(sce1.g, sce2.g,'stable');
        X1_ori = sce1.X(ia, :);
        X2_ori = sce2.X(ib, :);
    else
        g_ori = sce1.g;
        X1_ori = sce1.X;
        X2_ori = sce2.X;
    end
    
    X1_ori = sc_norm(X1_ori,'type','libsize');
    X2_ori = sc_norm(X2_ori,'type','libsize');

    [T1, X1, g1, xyz1] = sc_splinefit(X1_ori, g_ori, true, false);
    [T1, idx1] = sortrows(T1,'genes','ascend');
    X1 = X1(idx1, :);
    g1 = g1(idx1);

    [T2, X2, g2, xyz2] = sc_splinefit(X2_ori, g_ori, true, false);
    [T2, idx2] = sortrows(T2,'genes','ascend');
    X2 = X2(idx2, :);
    g2 = g2(idx2);

    assert(isequal(g1, g2))
    g = g1;

    px1 = T1.lgu; py1 = T1.lgcv; pz1 = T1.dropr;
    px2 = T2.lgu; py2 = T2.lgcv; pz2 = T2.dropr;

    %assignin("base","V1",[px1 py1 pz1]);
    %assignin("base","T1",T1);
    %assignin("base","xyz1",xyz1);

    v1=([px1 py1 pz1] - xyz1(T1.nearidx,:));
    v2=([px2 py2 pz2] - xyz2(T2.nearidx,:));

    DiffDist = vecnorm(v1 - v2, 2, 2);
    DiffSign = sign(vecnorm(v2,2,2)-vecnorm(v1,2,2));

    T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s', cL1{1}));
    T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s', cL2{1}));
    
    T = [T1 T2 table(DiffDist) table(DiffSign)];

    idxx = T.(8)==1 | T.(16)==1 | T.(8) == max(T.(8)) | T.(16) == max(T.(16));
    T.DiffDist(idxx) = 0;
    T = sortrows(T,"DiffDist","descend");
%}
    [T] = gui.e_dvanalysis(sce1, sce2, cL1, cL2);
    outfile = sprintf('%s_%s_vs_%s_%s.xlsx', ...
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
            'Sign of difference';...
            sprintf('%d',sce1.NumCells); sprintf('%d',sce2.NumCells)};
        Tnt = table(Item, Description);
        try
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end

        if strcmp(runenrichr,'Yes') && isok
            try
                [Tbp1, Tmf1]= run.r_enrichR(Tup.gene(1:min([250 size(Tup, 1)])));                
                in_writetable(Tbp1, filesaved, 'Up_250_GO_BP');
                in_writetable(Tmf1, filesaved, 'Up_250_GO_MF');
                [Tbp2, Tmf2]= run.r_enrichR(Tdn.gene(1:min([250 size(Tdn, 1)])));                
                in_writetable(Tbp2, filesaved, 'Dn_250_GO_BP');
                in_writetable(Tmf2, filesaved, 'Dn_250_GO_MF');
            catch ME
                disp(ME.message);
            end
        end
end
gui.gui_waitbar_adv(fw);

answer=questdlg(sprintf('Result files saved. Open the folder %s?', outdir), '');
if strcmp(answer,'Yes'), winopen(outdir); end

    function in_writetable(Tmf1, filesaved, shtname)
        if ~isempty(Tmf1) && istable(Tmf1) && size(Tmf1, 1) > 0
            writetable(Tmf1, filesaved, "FileType", "spreadsheet", 'Sheet', shtname);
        end
    end

end