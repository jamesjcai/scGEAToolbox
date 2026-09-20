function callback_GlycoShield(src, ~)
%GUI.CALLBACK_GLYCOSHIELD  Menu callback: cell-type-specific N-glycan
%shielding of docking targets.
%
%   Wraps GLY.SHIELD. Docking is run against PDB structures, which are
%   almost always solved deglycosylated or have their glycans stripped
%   during preparation. A compound can therefore score well against a bare
%   structure while never reaching the receptor on a real cell whose
%   glycocalyx covers the site.
%
%   The score is the product of two independent pieces of evidence: how
%   glycosylatable the protein is (N-X-S/T sequons in its UniProt sequence)
%   and how much N-glycosylation machinery the cell type actually runs
%   (normalized module scores). Neither alone is informative - a
%   sequon-rich protein in a cell that does not glycosylate is not
%   shielded, and neither is a sequon-free protein in one that does.
%
%   Sequences are fetched from UniProt, so this needs a network connection
%   the first time a gene is seen in a session.
%
% See also GLY.SHIELD, GLY.NORM, SC_DOCK_VINA, SC_DOCK_CCC, RUN.ML_SCDOCK.

[FigureHandle, sce] = gui.gui_getfigsce(src);

[labels, ~] = gui.i_getcellgroups(sce, FigureHandle);
if isempty(labels), return; end

[genes, affinity] = in_sourcegenes(FigureHandle);
if isempty(genes), return; end

% Restricting is worth offering rather than assuming: the output is one row
% per receptor per cell type, so all cell types on a long gene list is a
% large table for what is usually a question about one or two populations.
cellTypes = in_restrictcelltypes(labels, FigureHandle);
if isempty(cellTypes), return; end

nvargs = {};
if ~isempty(affinity)
    alpha = in_askalpha(FigureHandle);
    if isempty(alpha), return; end
    nvargs = {'affinity', affinity, 'alpha', alpha};
end
if ~isequal(sort(cellTypes), sort(unique(labels(labels ~= ""))))
    nvargs = [nvargs, {'cellType', cellTypes}];
end

answer = gui.myQuestdlg(FigureHandle, [ ...
    sprintf('Look up %d protein sequence(s) from UniProt?', numel(genes))
    ""
    "The sequon count comes from the sequence, so genes not already " + ...
    "cached this session are fetched over the network. Continue?"]);
if ~strcmp(answer, 'Yes'), return; end

fw = gui.myWaitbar(FigureHandle);
try
    T = gly.shield(genes, sce.X, sce.g, labels, nvargs{:});
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'gly.shield');
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(T)
    gui.myWarndlg(FigureHandle, ['No receptor could be scored. None of the ', ...
        'genes resolved to a UniProt sequence, or none of the cell types ', ...
        'survived the filter.']);
    return;
end

% A gene whose sequence never arrived is KEPT, with NaN in every
% sequence-derived column, rather than dropped. Those rows sort to the
% bottom and read like any other row at a glance, so name the genes: a
% missing UniProt lookup is a gap in the evidence, not a low score.
noSeq = unique(string(T.receptor(isnan(T.shield_score))));
msg = [sprintf('%d receptor x cell-type row(s).', height(T))
       ""
       "Shielding is a whole-protein estimate, not a pocket-level one: " + ...
       "SC_DOCK_VINA docks blind over a bounding box, so there is no " + ...
       "defined pocket to measure sequon distance from."];
if ~isempty(noSeq)
    msg(end+1) = "";
    msg(end+1) = sprintf(['%d gene(s) returned no UniProt sequence, so their ', ...
        'rows carry NaN rather than a shielding estimate: %s'], ...
        numel(noSeq), strjoin(in_head(noSeq, 6), ', '));
end
gui.myHelpdlg(FigureHandle, msg, 'N-glycan shielding');

gui.i_viewtable(T, FigureHandle);
gui.i_exporttable(T, true, 'Tglycoshield', 'GlycoShieldTable', ...
    [], [], FigureHandle);
end


function [genes, affinity] = in_sourcegenes(parentfig)
% Genes alone, or genes with the docking affinities to discount. Taking
% both from one table is what keeps them aligned: GLY.SHIELD wants the
% affinity vector element-for-element with the gene list, and pairing a
% pasted list against a separately chosen column is where that breaks.

genes = strings(0, 1);
affinity = [];

answer = gui.myQuestdlg(parentfig, [ ...
    "Which receptors should be scored?"
    ""
    "Paste a list: shielding scores only."
    "Workspace table: shielding scores plus discounted docking affinities, " + ...
    "taken from one table so the two stay aligned."], ...
    'Receptor genes', {'Paste a list', 'Workspace table', 'Cancel'}, 'Paste a list');

switch answer
    case 'Paste a list'
        glist = gui.i_inputgenelist([], false, parentfig);
        if isempty(glist), return; end
        genes = unique(string(glist), 'stable');
        genes = genes(strlength(genes) > 0);
    case 'Workspace table'
        T = gui.i_pickworkspacetable(parentfig, [], ...
            'Select a table with a receptor gene column and a docking affinity column:');
        if isempty(T), return; end
        [genes, affinity] = in_pickgenecolumns(T, parentfig);
    otherwise
        return;
end
end


function [genes, affinity] = in_pickgenecolumns(T, parentfig)
% One text column for the gene and one numeric column for the affinity.

genes = strings(0, 1);
affinity = [];
vars = string(T.Properties.VariableNames);

isnum = varfun(@isnumeric, T, OutputFormat = 'uniform');
textVars = vars(~isnum);
numVars = vars(isnum);
if isempty(textVars) || isempty(numVars)
    gui.myErrordlg(parentfig, ['The table needs both a text column of gene ', ...
        'symbols and a numeric column of affinities.'], 'gly.shield');
    return;
end

[indx, tf] = gui.myListdlg(parentfig, cellstr(textVars), ...
    'Which column holds the receptor gene symbols?', [], false);
if tf ~= 1 || isempty(indx), return; end
g = string(T.(textVars(indx)));

[indx, tf] = gui.myListdlg(parentfig, cellstr(numVars), ...
    'Which column holds the docking affinity (kcal/mol)?', [], false);
if tf ~= 1 || isempty(indx), return; end
a = double(T.(numVars(indx)));

keep = strlength(g) > 0 & ~ismissing(g) & isfinite(a);
% One row per gene: a duplicate would give the same receptor two different
% adjusted affinities in the output, with nothing to say which is meant.
[genes, ia] = unique(g(keep), 'stable');
aKeep = a(keep);
affinity = aKeep(ia);

if numel(genes) < sum(keep)
    gui.myWarndlg(parentfig, sprintf( ...
        ['%d duplicate gene row(s) dropped; the first affinity for each ', ...
         'gene was kept.'], sum(keep) - numel(genes)));
end
end


function sel = in_restrictcelltypes(labels, parentfig)
% Empty means cancelled; the full list means no restriction is passed on.

sel = strings(0, 1);
known = unique(labels(labels ~= ""));

answer = gui.myQuestdlg(parentfig, sprintf( ...
    ['Score all %d cell types?\n\nThe result is one row per receptor per ', ...
     'cell type, so restricting keeps the table readable.'], numel(known)));
switch answer
    case 'Yes'
        sel = known;
    case 'No'
        [indx, tf] = gui.myListdlg(parentfig, cellstr(known), ...
            'Select the cell types to score:', [], true);
        if tf ~= 1 || isempty(indx), return; end
        sel = known(indx);
    otherwise
        return;
end
end


function alpha = in_askalpha(parentfig)
% affinity_adj = affinity * (1 - alpha*shield_score). Vina affinities are
% negative and more negative is tighter, so shielding moves the value
% toward zero.

alpha = [];
answer = gui.i_inputdlg(['Strength of the affinity discount, 0 to 1 ', ...
    '(0 leaves affinities untouched, 1 applies the full shielding score):'], ...
    '0.5', parentfig);
if isempty(answer), return; end

v = str2double(answer{1});
if ~isfinite(v) || v < 0 || v > 1
    gui.myErrordlg(parentfig, 'Enter a number between 0 and 1.', 'gly.shield');
    return;
end
alpha = v;
end


function s = in_head(v, n)
% First N entries, with a count of what was left out.

v = string(v);
if numel(v) <= n
    s = v;
    return;
end
s = [v(1:n); sprintf('... and %d more', numel(v) - n)];
end
