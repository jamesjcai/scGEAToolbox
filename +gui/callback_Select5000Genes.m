function [requirerefresh] = callback_Select5000Genes(src)

requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

oldcn = sce.NumCells;
oldgn = sce.NumGenes;

if sce.NumGenes<=500
    gui.myWarndlg(FigureHandle, 'Number of genes is too small.');
    return;
end

spciestag = gui.i_selectspecies(2, false, FigureHandle);
if isempty(spciestag), return; end


prompt = {'Remove Mt-Genes (MT-ND1, MT-ND6, MT-CYB, MT-COI, MT-ATP6, etc.)?', ...
'Remove Hemoglobin Genes (HBA1, HBB, Hba-a1, etc.)?', ...
'Remove Genes With Name Contains ''orf'', ''-AS'', or ''-DT'' (C22orf42, C21orf58, etc.)?', ...
'Remove Genes With Name Starts With ''LINC'' (LINC01426, LINC01694, etc.)?', ...
'Remove Ribosomal Genes (RPSA, RPS2, RPS3, RPL3, RPL4, RPLP1, etc.)?', ...
'Remove Genes With Name Starts With ''Gm'' (Gm12768, Gm13305, etc.)?',...
'Remove Genes With Name Ends With ''Rik'' (0610005C13Rik, 0610007C21Ri, etc.)?',...
'Remove Genes Without Approved Symbols?', ...
'Remove Genes Expressed in Less Than m Cells (m = 0.075 or 0.050, 10 or 50)?', ...
'Keep Top n Highly Variable Genes (HVGs) (n = 5000 or 2000)?', ...
'Whitelist - Genes to Keep No Matter What (Cblb, Lyz2, ...; Blank For None)?'};
dlgtitle = '';
dims = [1, 80];

definput = {'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', '0.075', num2str(min([sce.NumGenes,5000])), ''};

% Field positions, named rather than counted from the end. The two numeric
% answers used to be read as answer{end-1} and answer{end}; appending the
% whitelist field would have silently re-pointed those at the wrong boxes.
IDX_MINCELLS  = 9;
IDX_HVGN      = 10;
IDX_WHITELIST = 11;


if gui.i_isuifig(FigureHandle)
    answer = gui.myInputdlg(prompt, dlgtitle, definput, FigureHandle);
else
    answer = inputdlg(prompt, dlgtitle, dims, definput);
end
if isempty(answer)
    requirerefresh = false;
    return;
end


% ---- whitelist ----------------------------------------------------------
% Genes the user never wants dropped, whatever the filters decide. Resolved
% against sce.g now, while every gene is still present, and stashed with
% their expression so anything a later step removes can be put back. This is
% the same stash-and-restore idiom SingleCellExperiment.qcfilterwhitelist
% uses internally.
[whitelist, unknown] = i_parsewhitelist(answer{IDX_WHITELIST}, sce.g);
if ~isempty(unknown)
    gui.myWarndlg(FigureHandle, sprintf(['These whitelist genes are not in ' ...
        'this dataset and will be ignored:\n%s'], strjoin(unknown, ', ')));
end
Xwhite = [];
if ~isempty(whitelist)
    [~, wi] = ismember(whitelist, sce.g);
    Xwhite = sce.X(wi, :);
    fprintf('%d whitelisted gene(s) will be kept regardless of filters.\n', ...
        numel(whitelist));
end

if sce.NumCells*sce.NumGenes < 4e8
    sceori = copy(sce);
    % disp('Ready for reversible.');
else
    confirmanswer = gui.myQuestdlg(FigureHandle, 'You are about to change the SCE data. This cannot be undone.');
    if ~strcmp(confirmanswer, 'Yes'), return; end
    sceori = [];
end

fw = gui.myWaitbar(FigureHandle);

c = 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    sce = sce.rmmtgenes;
    disp('Mt-genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    sce = sce.rmhemoglobingenes;
    disp('Hemoglobin genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = contains(sce.g, 'orf') | contains(sce.g, '-AS') | contains(sce.g, '-as') | contains(sce.g, '-DT') | contains(sce.g, '-dt');
    sce.X(idx, :) = [];
    sce.g(idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes with name contains ''orf'' or ''-AS'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = startsWith(sce.g, 'LINC');
    sce.X(idx, :) = [];
    sce.g(idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes with name starts with ''LINC'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    sce = sce.rmribosomalgenes;
    disp('Ribosomal genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = find(~cellfun(@isempty, regexp(sce.g,"Gm[0-9][0-9][0-9]")));
    sce.X(idx, :) = [];
    sce.g(idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes with name starts with ''Gm'' are found and removed.\n',a1-a2);
end


c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = endsWith(sce.g, 'Rik');
    sce.X(idx, :) = [];
    sce.g(idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes with name ends with ''Rik'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    mfolder = fileparts(mfilename('fullpath'));
    switch spciestag
        case 'human'
            % T = readtable(fullfile(mfolder, ...
            %     '../assets', 'HGNCBiomart.txt'));
            load(fullfile(mfolder, ...
                '..', 'assets', 'Biomart', 'Biomart_human_genes.mat'), 'T');
        case 'mouse'
            load(fullfile(mfolder, ...
                '..', 'assets',  'Biomart', 'Biomart_mouse_genes.mat'), 'T');
    end
    ApprovedSymbol = string(T.GeneName);
    [idx] = ismember(upper(sce.g), upper(ApprovedSymbol));
    a1 = length(sce.g);
    sce.X(~idx, :) = [];
    sce.g(~idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes without approved symbols are found and removed.\n',a1-a2);
    % requirerefresh = true;
end

try
    a = str2double(answer{IDX_MINCELLS});
    if a > 0 && a < intmax
        a1 = length(sce.g);
        sce = sce.selectkeepgenes(1, a);
        a2 = length(sce.g);
        fprintf('%d lowly expressed genes found and removed.\n',a1-a2);
        % requirerefresh = true;
    end
catch ME
    warning(ME.message);
end

try
    a = str2double(answer{IDX_HVGN});
    % T = sc_hvg(sce.X, sce.g);
    T = sc_splinefit(sce.X, sce.g);
    glist = T.genes(1:min([a, sce.NumGenes]));
    [y, idx] = ismember(glist, sce.g);
    if ~all(y)
        gui.myErrordlg(FigureHandle, 'Runtime error.');
        return;
    end
    sce.X = sce.X(idx, :);
    sce.g = sce.g(idx);
 catch ME
     gui.myWaitbar(FigureHandle, fw,true);
     gui.myWarndlg(FigureHandle, ME.message, ME.identifier);
     return;
 end

% ---- put back anything the gene filters dropped --------------------------
% Every step above removes genes only, never cells, so the stashed rows
% still line up column-for-column and can simply be appended. Doing this
% before the QC call also satisfies that method's assert, which requires
% each whitelisted gene to be present in sce.g when it is handed over.
if ~isempty(whitelist)
    missing = ~ismember(whitelist, sce.g);
    if any(missing)
        sce.X = [sce.X; Xwhite(missing, :)];
        sce.g = [sce.g; whitelist(missing)];
        fprintf('%d whitelisted gene(s) restored after filtering: %s\n', ...
            sum(missing), strjoin(whitelist(missing), ', '));
    end
end

try
    % assignin("base", "c_before", sce.c);
    % assignin("base", "numcells_before", sce.numcells);
    % The whitelist goes through to the QC step too - that one removes
    % CELLS, so it does its own stash-and-restore with the surviving
    % columns rather than the ones captured up top.
    sce = sce.qcfilterwhitelist(1000, 0.15, 15, 500, whitelist);
    % assignin("base", "c_after", sce.c);
    % assignin("base", "numcells_after", sce.numcells);
catch ME
    gui.myWaitbar(FigureHandle, fw,true);
    gui.myWarndlg(FigureHandle, ME.message, ME.identifier);
    return;
end

gui.myWaitbar(FigureHandle, fw);

newcn = sce.NumCells;
newgn = sce.NumGenes;

% sce is a handle object that was modified in place, so restoring sceori
% requires pushing the pristine copy back to the caller.
if newgn==0
    if ~isempty(sceori)
        gui.myHelpdlg(FigureHandle, "All genes are removed. Operation is cancelled.");
        sce = copy(sceori);
        gui.myGuidata(FigureHandle, sce, src);
    else
        requirerefresh = true;
    end
    return;
end
if newcn==0
    if ~isempty(sceori)
        gui.myHelpdlg(FigureHandle, "All cells are removed. Operation is cancelled.");
        sce = copy(sceori);
        gui.myGuidata(FigureHandle, sce, src);
    else
        requirerefresh = true;
    end
    return;
end
if oldcn-newcn==0 && oldgn-newgn==0
    gui.myHelpdlg(FigureHandle, "No cells and genes are removed.");
    return;
end
if ~isempty(sceori)
    acceptanswer = gui.myQuestdlg(FigureHandle, sprintf('%d genes will be removed; %d cells will be removed.\n[%d genes x %d cells] => [%d genes x %d cells]', ...
            oldgn-newgn, oldcn-newcn, oldgn, oldcn, newgn, newcn),'', ...
            {'Accept Changes', 'Cancel Changes'}, 'Accept Changes');
    if ~strcmp(acceptanswer, 'Accept Changes')
        sce = copy(sceori);
    else
        requirerefresh = true;
    end
else
    requirerefresh = true;
end
gui.myGuidata(FigureHandle, sce, src);

end


function [whitelist, unknown] = i_parsewhitelist(txt, glist)
%I_PARSEWHITELIST  Free-text gene names -> names as spelled in glist.
%
% Accepts commas, semicolons, whitespace or newlines as separators, so a
% list pasted from a spreadsheet, a paper or a previous MATLAB session all
% work without reformatting. Matching is case-insensitive but the returned
% names use the dataset's own spelling, because every downstream ismember
% against sce.g is case-sensitive.
whitelist = strings(0, 1);
unknown   = strings(0, 1);
if isempty(txt), return; end

parts = strtrim(split(string(txt), [",", ";", " ", newline, char(9)]));
parts = parts(strlength(parts) > 0);
if isempty(parts), return; end
parts = unique(parts, 'stable');

g = string(glist(:));
[tf, loc] = ismember(upper(parts), upper(g));
whitelist = g(loc(tf));                 % dataset spelling, not the user's
whitelist = unique(whitelist, 'stable');
unknown   = parts(~tf);
end
