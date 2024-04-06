function callback_DEGene2GroupsBatch(src, ~)


FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);


if length(unique(sce.c_cell_type_tx))==1
    warndlg('Only one cell type or cell type is undetermined.','');
    return;
end

if ~gui.gui_showrefinfo('DE in Batch Mode'), return; end

[CellTypeSorted] = pkg.e_sortcatbysize(sce.c_cell_type_tx);
[CellTypeList] = in_selectcelltypes(CellTypeSorted);
if isempty(CellTypeList), return; end

[thisc, clabel]=in_select1class(sce);
if isempty(thisc), return; end
if strcmp(clabel,'Cell Type')
    helpdlg('Cannot select ''Cell Type'' as grouping varialbe.');
    return;
end

%[i1, i2, cL1, cL2] = gui.i_select2grps(sce, false);

[i1, i2, cL1, cL2]=in_twogrpsencoding(thisc);
if isempty(i1) || isempty(i2) || isempty(cL1) || isempty(cL2)
    return;
end
if length(i1) == 1 || length(i2) == 1, return; end

%methodtag = "ranksum";

% answer = questdlg('Which method?', ...
%     'Select Method', 'Wilcoxon rank-sum test 🐇', ...
%     'DESeq 2 (R required) 🐢', ...
%     'MAST (R required) 🐢', ...
%     'Wilcoxon rank-sum test 🐇');
% 
% if strcmpi(answer, 'Wilcoxon rank-sum test 🐇')
%     methodtag = "ranksum";
% elseif strcmpi(answer, 'DESeq 2 (R required) 🐢')
%     methodtag = "deseq2";
% elseif strcmpi(answer, 'MAST (R required) 🐢')
%     methodtag = "mast";
%     if isempty(pkg.FindRpath)
%         warndlg('This function requires R environment.')
%             return;
%     end
% else
%     return;
% end
answer=questdlg('Select a folder to save the outupt Excel files. Continue?','');
if ~strcmp(answer,'Yes'), return; end

outdir = uigetdir;
if ~isfolder(outdir), return; end

needoverwritten=false;
for k=1:length(CellTypeList)
    outfile = sprintf('%s_vs_%s_%s.xlsx', ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)), ...
        matlab.lang.makeValidName(string(CellTypeList{k})));
    filesaved = fullfile(outdir, outfile);
    if exist(filesaved,'file')
        needoverwritten=true;
    end
end
if needoverwritten
    answer=questdlg(sprintf('Overwrite existing result file(s) in %s?',outdir),'');    
else
    answer=questdlg(sprintf('Result files will be save in %s. Continue?', outdir), '');
end
if ~strcmp(answer,'Yes'), return; end

rungsea=false;
answer=questdlg('Run GSEA analysis with the DE result (i.e., genes ranked by log2FC)?','');
if strcmp(answer,'Cancel'), return; end
if strcmp(answer,'Yes'), rungsea = true; end

fw = gui.gui_waitbar_adv;

for k=1:length(CellTypeList)
   
    gui.gui_waitbar_adv(fw, ...
        (k-1)/length(CellTypeList), ...
        sprintf('Processing %s ...', CellTypeList{k}));

    outfile = sprintf('%s_vs_%s_%s.xlsx', ...
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

    
    T = sc_deg(sce.X(:, i1&idx), sce.X(:, i2&idx), sce.g, 1, false);
    T = in_DETableProcess(T,cL1,cL2);
    [Tup, Tdn] = pkg.e_processDETable(T);
    try
        writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
        writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
        writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
    catch ME
        warning(ME.message);
    end
    if rungsea
        try
            Tgsea=run.r_fgsea(T.gene,true,T.avg_log2FC);
            if ~isempty(Tgsea)
                writetable(Tgsea, filesaved, "FileType", "spreadsheet", 'Sheet', 'GSEA');
            end
        catch ME
            warning(ME.message);
        end
    end
end
gui.gui_waitbar_adv(fw);

answer=questdlg(sprintf('Result files saved. Open the folder %s?', outdir), '');
if strcmp(answer,'Yes'), winopen(outdir); end

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






function [i1, i2, cL1, cL2]=in_twogrpsencoding(thisc)
    [ci, cLi] = grp2idx(thisc);
    listitems = natsort(string(cLi));
    n = length(listitems);
    if n < 2
        errordlg('Need at least two groups.');
        return;
    end
    [indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listitems, ...
        'InitialValue', [n - 1, n], ...
        'ListSize', [220, 300]);
    if tfx == 1
        if numel(indxx) ~= 2
            errordlg('Please select 2 groups');
            return;
        end
        [y1, idx1] = ismember(listitems(indxx(1)), cLi);
        [y2, idx2] = ismember(listitems(indxx(2)), cLi);
        assert(y1 & y2);
        i1 = ci == idx1;
        i2 = ci == idx2;
        cL1 = cLi(idx1);
        cL2 = cLi(idx2);
    else
        i1=[]; i2=[]; cL1=[]; cL2=[];
    end
end

function [CellTypeList]=in_selectcelltypes(CellTypeSorted)
    CellTypeList=[];
    %pause(1);
    %[idx] = gui.i_selmultidlg(CellTypeSorted, CellTypeSorted);
    %if isempty(idx), return; end
    %if idx == 0, return; end
    %CellTypeList = CellTypeSorted(idx);
    [indx2, tf2] = listdlg('PromptString', ...
    {'Select Cell Types:'}, ...
    'SelectionMode', 'multiple', 'ListString', CellTypeSorted, ...
    'InitialValue',1:length(CellTypeSorted) ...
    'ListSize', [220, 300]);
    if tf2 == 1
        CellTypeList = CellTypeSorted(indx2);
    end
end

% function [cLsorted]=in_sortcelltypes(thisc)
%     [c, cL] = grp2idx(thisc);
%     cmv = 1:max(c);
%     %idxx = cmv;
%     [cmx] = countmember(cmv, c);
%     [~, idxx] = sort(cmx, 'descend');
%     cLsorted=cL(idxx);
% end

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



    function [thisc, clabel] = in_select1class(sce, allowunique)
        if nargin < 2, allowunique = false; end
        thisc = [];
        clabel = '';
        
        listitems = {'Current Class (C)'};
        if ~isempty(sce.c_cluster_id)
            if allowunique
                listitems = [listitems, 'Cluster ID'];
            else
                if numel(unique(sce.c_cluster_id)) > 1
                    listitems = [listitems, 'Cluster ID'];
                end
            end
        end
        
        if ~isempty(sce.c_cell_cycle_tx)
            if allowunique
                listitems = [listitems, 'Cell Cycle Phase'];
            else
                if numel(unique(sce.c_cell_cycle_tx)) > 1
                    listitems = [listitems, 'Cell Cycle Phase'];
                end
            end
        end
        if ~isempty(sce.c_batch_id)
            if allowunique
                listitems = [listitems, 'Batch ID'];
            else
                if numel(unique(sce.c_batch_id)) > 1
                    listitems = [listitems, 'Batch ID'];
                end
            end
        end
        
        a = evalin('base', 'whos');
        b = struct2cell(a);
        v = false(length(a), 1);
        for k = 1:length(a)
            if max(a(k).size) == sce.NumCells && min(a(k).size) == 1
                v(k) = true;
            end
        end
        if any(v)
            a = a(v);
            b = b(:, v);
            listitems = [listitems, 'Workspace Variable...'];
        end
        
        [indx2, tf2] = listdlg('PromptString', ...
            {'Select grouping variable:'}, ...
            'SelectionMode', 'single', ...
            'ListString', listitems ...
            'ListSize', [220, 300]);
        if tf2 == 1
            clabel = listitems{indx2};
            switch clabel
                case 'Current Class (C)'
                    thisc = sce.c;
                case 'Cluster ID' % cluster id
                    thisc = sce.c_cluster_id;
                case 'Batch ID' % batch id
                    thisc = sce.c_batch_id;
                %case 'Cell Type' % cell type
                %    thisc = sce.c_cell_type_tx;
                case 'Cell Cycle Phase' % cell cycle
                    thisc = sce.c_cell_cycle_tx;
                case 'Workspace Variable...'
                    thisc = i_pickvariable;
            end
        end
        
        function [c] = i_pickvariable
            c = [];
            [indx, tf] = listdlg('PromptString', {'Select variable:'}, ...
                'liststring', b(1, :), ...
                'SelectionMode', 'single', ...
                'ListSize', [220, 300]);
            if tf == 1
                c = evalin('base', a(indx).name);
            end
            %    end
        end
    end
