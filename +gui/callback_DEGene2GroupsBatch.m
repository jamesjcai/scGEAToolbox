function callback_DEGene2GroupsBatch(src, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
if ~gui.gui_showrefinfo('DE in Batch Mode'), return; end
prefixtag = 'DE';

[done, CellTypeList, i1, i2, cL1, cL2,... 
    outdir] = gui.i_batchmodeprep(sce, prefixtag);
if ~done, return; end

runenrichr = questdlg('Run Enrichr (R required) with top 250 DE genes? Results will be saved in the output Excel files.','');
if strcmp(runenrichr,'Cancel'), return; end

[isok, msg] = commoncheck_R('R_enrichR');
if ~isok
    disp(msg);
    answer = questdlg('R Environment Error: It seems that your R environment is not configured correctly. This means that Gene Function Enrichment Analysis using enrichR cannot be performed for differentially expressed genes. Continue withouth enrichR?',''); 
    if ~strcmp(answer,'Yes'), return; end
end

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
    % if exist(filesaved,'file')
    %     answer=questdlg(sprintf('Overwrite existing file %s? Click No to skip.',outfile),'');
    %     switch answer
    %         case 'Yes'
    %         case 'No'
    %             continue;
    %         case 'Cancel'
    %             return;
    %     end
    % end

    idx = sce.c_cell_type_tx == CellTypeList{k};
    T = [];
    try
        T = sc_deg(sce.X(:, i1&idx), sce.X(:, i2&idx), sce.g, 1, false);
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
            'adjusted p-value';...
             sprintf('%d',sum(i1&idx)); sprintf('%d',sum(i2&idx))};
        Tnt = table(Item, Description);
%{
gene
p_val
avg_log2FC
abs_log2FC
avg_1_GSM3308547
avg_2_GSM3308548
pct_1_GSM3308547
pct_2_GSM3308548
p_val_adj
%}
        
        [Tup, Tdn] = pkg.e_processDETable(T);
        try
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
            writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
        catch ME
            warning(ME.message);
        end

        if strcmp(runenrichr,'Yes')
            try
                [Tbp1, Tmf1]= run.r_enrichR(Tup.gene(1:min([250 size(Tup, 1)])));
                %[Tbp1, Tmf1]= run.py_GSEApy_enr(Tup.gene(1:min([250 size(Tup, 1)])), T.gene, tempdir);
                in_writetable(Tbp1, filesaved, 'Up_250_GO_BP');
                in_writetable(Tmf1, filesaved, 'Up_250_GO_MF');                
                [Tbp2, Tmf2]= run.r_enrichR(Tdn.gene(1:min([250 size(Tdn, 1)])));
                %[Tbp2, Tmf2]= run.py_GSEApy_enr(Tdn.gene(1:min([250 size(Tup, 1)])), T.gene, tempdir);                
                in_writetable(Tbp2, filesaved, 'Dn_250_GO_BP');
                in_writetable(Tmf2, filesaved, 'Dn_250_GO_MF');
            catch ME
                disp(ME.message);
            end
        end
        
        % if rungsea
        %     try
        %         Tgsea=run.r_fgsea(T.gene,true,T.avg_log2FC);
        %         if ~isempty(Tgsea)
        %             writetable(Tgsea, filesaved, "FileType", "spreadsheet", 'Sheet', 'GSEA');
        %         end
        %     catch ME
        %         warning(ME.message);
        %     end
        % end
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


function [T]=in_DETableProcess(T,cL1,cL2)
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


% try
%     switch methodtag
%         case 'ranksum'
%             %fw=gui.gui_waitbar;
%             T = sc_deg(sce.X(:, i1), sce.X(:, i2), ...
%                 sce.g, 1, true);
%         case 'deseq2'
%             [ok] = gui.i_confirmscript('DE analysis (DESeq2)', ...
%                 'R_DESeq2', 'r');
%             if ~ok, return; end
%             fw = gui.gui_waitbar;
%             T = run.r_DESeq2(sce.X(:, i1), sce.X(:, i2), sce.g);
%             gui.gui_waitbar(fw);
%         case 'mast'
%             [ok] = gui.i_confirmscript('DE analysis (MAST)', ...
%                 'R_MAST', 'r');
%             if ~ok, return; end
%             fw = gui.gui_waitbar;
%             T = run.r_MAST(sce.X(:, i1), sce.X(:, i2), sce.g);
%             gui.gui_waitbar(fw);
%     end
% 
% catch ME
%     gui.gui_waitbar(fw, true);
%     errordlg(ME.message);
%     return;
% end

% mavolcanoplot(sce.X(:,i1),sce.X(:,i2),T.p_val_adj,'Labels',T.gene)
