function callback_DEGene2GroupsBatch(src, ~)


FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

warndlg('Function is under development.','');
return;


answer = questdlg('Compare 2 groups across all cell types (DE analysis in batch)?', '');
if ~strcmp(answer,'Yes'), return; end

[thisc, clable]=gui.i_select1class(sce);
if strcmp(clable,'Cell Type')
    helpdlg('Cannot select ''Cell Type'' as grouping varialbe.');
    return;
end

%[i1, i2, cL1, cL2] = gui.i_select2grps(sce);

[i1, i2, cL1, cL2]=in_twogrpsencoding(thisc);
if isempty(i1) || isempty(i2) || isempty(cL1) || isempty(cL2)
    return;
end
if length(i1) == 1 || length(i2) == 1, return; end

answer = questdlg('Which method?', ...
    'Select Method', 'Wilcoxon rank-sum test 🐇', ...
    'DESeq 2 (R required) 🐢', ...
    'MAST (R required) 🐢', ...
    'Wilcoxon rank-sum test 🐇');

if strcmpi(answer, 'Wilcoxon rank-sum test 🐇')
    methodtag = "ranksum";
elseif strcmpi(answer, 'DESeq 2 (R required) 🐢')
    methodtag = "deseq2";
elseif strcmpi(answer, 'MAST (R required) 🐢')
    methodtag = "mast";
    if isempty(pkg.FindRpath)
        warndlg('This function requires R environment.')
            return;
    end
else
    return;
end

% ------------
[CellTypeSorted]=in_sortcelltypes(sce.c_cell_type_tx)

for k=1:length(CellTypeSorted)
    idx1 = sce.c_cell_type_tx == CellTypeSorted{k};
    % sum(idx1)
end



try
    switch methodtag
        case 'ranksum'
            %fw=gui.gui_waitbar;
            T = sc_deg(sce.X(:, i1), sce.X(:, i2), ...
                sce.g, 1, true);
        case 'deseq2'
            [ok] = gui.i_confirmscript('DE analysis (DESeq2)', ...
                'R_DESeq2', 'r');
            if ~ok, return; end
            fw = gui.gui_waitbar;
            T = run.r_DESeq2(sce.X(:, i1), sce.X(:, i2), sce.g);
            gui.gui_waitbar(fw);
        case 'mast'
            [ok] = gui.i_confirmscript('DE analysis (MAST)', ...
                'R_MAST', 'r');
            if ~ok, return; end
            fw = gui.gui_waitbar;
            T = run.r_MAST(sce.X(:, i1), sce.X(:, i2), sce.g);
            gui.gui_waitbar(fw);
    end

catch ME
    gui.gui_waitbar(fw, true);
    errordlg(ME.message);
    return;
end

% mavolcanoplot(sce.X(:,i1),sce.X(:,i2),T.p_val_adj,'Labels',T.gene)

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

outfile = sprintf('%s_vs_%s', ...
    matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
if isatac, T.gene = "chr" + T.gene; end

[filetype, filesaved] = gui.i_exporttable(T, true, 'T', outfile);
tf = 0;
if ~(ismcc || isdeployed) && strcmp(filetype, 'Workspace')
    [Tup, Tdn] = pkg.e_processDETable(T, true);
    labels = {'Save DE results (selected up-regulated) to variable named:', ...
        'Save DE results (selected down-regulated) to variable named:'};
    vars = {'Tup', 'Tdn'};
    values = {Tup, Tdn};
    [~, tf] = export2wsdlg(labels, vars, values);
    if tf ~= 1
        return;
    end
end

if ~isempty(filesaved)
    if strcmp(filetype, 'Excel file')
        answer = questdlg('Save up- and down-regulated genes to seperate sheets?');
        if strcmp(answer, 'Yes')
            [Tup, Tdn] = pkg.e_processDETable(T, true);
            % strcmp(extractAfter(filesaved,strlength(filesaved)-4),'xlsx')
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
            waitfor(helpdlg(sprintf('Result has been saved in %s', filesaved), ''));
            %writetable(Tup,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet',);
            %writetable(Tdn,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet');
        end
    elseif strcmp(filetype, 'Text file')
        % strcmp(extractAfter(filesaved,strlength(filesaved)-3),'txt')
        answer = questdlg('Save up- and down-regulated genes to seperate files?');
        if strcmp(answer, 'Yes')
            [Tup, Tdn] = pkg.e_processDETable(T, true);
            [~, ~] = gui.i_exporttable(Tup, true, 'Tup', 'Upregulated', 'Text file');
            [~, ~] = gui.i_exporttable(Tdn, true, 'Tdn', 'Downregulated', 'Text file');
        end
    end
end


    if tf == 1
        disp('To run pathway analysis, type:');
        disp('t=run.r_SPIA(T);');
        disp('To run enrichment analysis, type:');
        disp('run.web_Enrichr(Tup.gene(1:200))');
        disp('run.web_Enrichr(Tdn.gene(1:200))');

        answer = questdlg('Run enrichment analysis with top K (=200 by default) up-regulated DE genes?');
        if strcmp(answer, 'Yes')
            gui.i_enrichtest(Tup.gene(1:min(numel(Tup.gene), 200)), sce.g);
        else
            return;
        end

        answer = questdlg('Run enrichment analysis with top K (=200 by default) down-regulated DE genes?');
        if strcmp(answer, 'Yes')
            gui.i_enrichtest(Tdn.gene(1:min(numel(Tdn.gene), 200)), sce.g);
        else
            return;
        end
    end
end

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
        'InitialValue', [n - 1, n]);
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

function [cLsorted]=in_sortcelltypes(thisc)
    [c, cL] = grp2idx(thisc);
    cmv = 1:max(c);
    %idxx = cmv;
    [cmx] = countmember(cmv, c);
    [~, idxx] = sort(cmx, 'descend');
    cLsorted=cL(idxx);
end
