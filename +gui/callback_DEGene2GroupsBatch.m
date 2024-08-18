function callback_DEGene2GroupsBatch(src, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
if ~gui.gui_showrefinfo('DE in Batch Mode'), return; end
prefixtag = 'DE';

[done, CellTypeList, i1, i2, cL1, cL2,... 
    outdir] = gui.i_batchmodeprep(sce, prefixtag);
if ~done, return; end

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
        warning(ME.message);
    end

    if ~isempty(T)
        T = in_DETableProcess(T, cL1, cL2);
        [Tup, Tdn] = pkg.e_processDETable(T);
        try
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
        catch ME
            warning(ME.message);
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
