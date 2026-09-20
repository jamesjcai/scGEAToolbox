function callback_GlycoWeight(src, ~)
%GUI.CALLBACK_GLYCOWEIGHT  Menu callback: re-weight cell-cell communication
%edges by glyco-state.
%
%   Wraps GLY.WEIGHT, which discounts an existing PROTEIN ligand-receptor
%   table by the glycan context each interaction needs - heparan sulfate for
%   FGF/Wnt/chemokines, Fringe glycosylation for Notch, N-glycan branching
%   for EGFR-family residency. Edges with no glyco dependency keep a factor
%   of 1, and the original probability column is never overwritten.
%
%   The table has to come from somewhere else: SC_DOCK_CCC produces one, and
%   it has no menu entry, so this reads the table from the base workspace or
%   from a file. GUI.CALLBACK_GLYCOCCC does NOT produce a suitable table -
%   its rows are glycan-to-lectin pairs with no ligand/receptor columns.
%
% See also GLY.WEIGHT, GLY.LECTINMAP, GLY.NORM, SC_DOCK_CCC,
% GUI.CALLBACK_GLYCOCCC.

[FigureHandle, sce] = gui.gui_getfigsce(src);

T_ccc = in_sourcetable(FigureHandle);
if isempty(T_ccc), return; end

vars = string(T_ccc.Properties.VariableNames);
missing = setdiff(["ligand", "receptor", "sender", "receiver"], vars, 'stable');
if ~isempty(missing)
    gui.myErrordlg(FigureHandle, sprintf( ...
        ['The table is missing the column(s): %s.\n\n', ...
         'This analysis re-weights a protein ligand-receptor table - one ', ...
         'row per directed edge, with ligand, receptor, sender and ', ...
         'receiver columns, as SC_DOCK_CCC returns in T_interactions. A ', ...
         'glyco-lectin table from the Glyco-Lectin Cell-Cell Communication ', ...
         'entry has glyco_module and lectin instead and cannot be used here.'], ...
        strjoin(missing, ', ')), 'gly.weight');
    return;
end

[labels, ~] = gui.i_getcellgroups(sce, FigureHandle);
if isempty(labels), return; end

% A table built on a different run carries cell-type names this dataset does
% not have. Every edge then falls through to a factor of 1, and the output
% looks like "no glyco dependency anywhere" rather than like a mismatch.
edgeTypes = unique([string(T_ccc.sender); string(T_ccc.receiver)]);
known = unique(labels(labels ~= ""));
if ~any(ismember(edgeTypes, known))
    gui.myErrordlg(FigureHandle, sprintf( ...
        ['None of the cell types named in the table (%s) appear in this ', ...
         'dataset (%s).\n\nThe weights are looked up per cell type, so ', ...
         'every edge would silently get a factor of 1. Load the dataset ', ...
         'the table was built from, or pick the grouping whose names it ', ...
         'uses.'], ...
        strjoin(in_head(edgeTypes, 4), ', '), strjoin(in_head(known, 4), ', ')), ...
        'gly.weight');
    return;
elseif ~all(ismember(edgeTypes, known))
    unknown = edgeTypes(~ismember(edgeTypes, known));
    answer = gui.myQuestdlg(FigureHandle, sprintf( ...
        ['%d of %d cell types named in the table are not in this dataset ', ...
         '(%s).\n\nThose edges will keep a factor of 1 whatever their ', ...
         'glycan dependency. Continue?'], numel(unknown), numel(edgeTypes), ...
        strjoin(in_head(unknown, 4), ', ')));
    if ~strcmp(answer, 'Yes'), return; end
end

% A two-condition table's probability column is a between-condition
% DIFFERENCE. Weighting it with one condition-averaged factor cancels out
% exactly the glycan remodelling the comparison is for - GLY.WEIGHT warns
% about this on the console, which is nowhere a GUI user will look, so make
% it a decision instead.
if all(ismember(["prob_cond1", "prob_cond2"], vars))
    [cond, cond1, cond2] = in_askcondition(sce, FigureHandle);
    if isempty(cond), return; end
    nvargs = {'condition', cond, 'condLabels', [cond1; cond2]};
else
    probVar = in_askprobvar(T_ccc, vars, FigureHandle);
    if strlength(probVar) == 0, return; end
    nvargs = {'probVar', probVar};
end

[~, methodid] = gui.i_pickscoremethod([], FigureHandle);
if isempty(methodid), return; end

fw = gui.myWaitbar(FigureHandle);
try
    T_out = gly.weight(T_ccc, sce.X, sce.g, labels, nvargs{:}, methodid=methodid);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'gly.weight');
    return;
end
gui.myWaitbar(FigureHandle, fw);

nModulated = sum(strlength(string(T_out.glyco_module)) > 0);
gui.myHelpdlg(FigureHandle, [ ...
    sprintf('%d of %s re-weighted after matching a glyco dependency.', ...
        nModulated, pkg.i_plural(height(T_out), 'edge'))
    ""
    "The factor is min-max normalized across the groups present, so 0 is " + ...
    "the lowest-scoring group and 1 the highest. It is a relative " + ...
    "statement about this dataset, not an absolute amount of glycan."], ...
    'Glyco re-weighted communication');

gui.i_viewtable(T_out, FigureHandle);
gui.i_exporttable(T_out, true, 'Tglycoweight', 'GlycoWeightTable', ...
    [], [], FigureHandle);
end


function T = in_sourcetable(parentfig)
% The table cannot come from the SCE, so ask where it does come from.

T = [];
answer = gui.myQuestdlg(parentfig, [ ...
    "Where is the ligand-receptor communication table?"
    ""
    "It needs ligand, receptor, sender and receiver columns plus a " + ...
    "probability column, as SC_DOCK_CCC returns in T_interactions."], ...
    'Communication table', {'Base workspace', 'File', 'Cancel'}, 'Base workspace');
switch answer
    case 'Base workspace'
        T = gui.i_pickworkspacetable(parentfig, [], ...
            'Select the ligand-receptor communication table:');
    case 'File'
        T = gui.i_readtablefile(parentfig, [], ...
            'Pick a ligand-receptor communication table');
    otherwise
        return;
end
end


function probVar = in_askprobvar(T_ccc, vars, parentfig)
% "prob" is the default name and the usual one; when it is absent, the
% column has to be named rather than guessed at.

probVar = "";
if ismember("prob", vars)
    probVar = "prob";
    return;
end

numericVars = vars(varfun(@isnumeric, T_ccc, OutputFormat = 'uniform'));
if isempty(numericVars)
    gui.myErrordlg(parentfig, ['The table has no numeric column to use as ', ...
        'the edge probability.'], 'gly.weight');
    return;
end

[indx, tf] = gui.myListdlg(parentfig, cellstr(numericVars), ...
    'Which column holds the edge probability?', [], false);
if tf ~= 1 || isempty(indx), return; end
probVar = numericVars(indx);
end


function [cond, cond1, cond2] = in_askcondition(sce, parentfig)
% Which per-cell attribute the table's two conditions correspond to, and
% which of its values is prob_cond1. The order is not cosmetic: it decides
% which factor multiplies which probability.

cond = strings(0, 1);
cond1 = "";
cond2 = "";

gui.myHelpdlg(parentfig, [ ...
    "This table carries prob_cond1 and prob_cond2, so it came from a " + ...
    "two-condition run and its probability column is a between-condition " + ...
    "difference."
    ""
    "Glycosylation programs are themselves remodelled between conditions, " + ...
    "so each condition has to be weighted by its own factor. Applying one " + ...
    "averaged factor to a difference cancels out the remodelling the " + ...
    "comparison is meant to detect."
    ""
    "Select the per-cell attribute holding the two conditions."], ...
    'Two-condition table');

[thisc, clabel] = gui.i_select1class(sce, false, ...
    'Select the condition attribute:', '', parentfig);
if isempty(thisc), return; end

vals = unique(string(thisc), 'stable');
vals = vals(vals ~= "");
if numel(vals) ~= 2
    gui.myErrordlg(parentfig, sprintf( ...
        ['"%s" has %d distinct values, but a two-condition table needs ', ...
         'exactly 2.'], clabel, numel(vals)), 'gly.weight');
    return;
end

[indx, tf] = gui.myListdlg(parentfig, cellstr(vals), ...
    'Which condition is prob_cond1?', [], false);
if tf ~= 1 || isempty(indx), return; end

cond = string(thisc);
cond1 = vals(indx);
cond2 = vals(setdiff(1:2, indx));
end


function s = in_head(v, n)
% First N entries, with a count of what was left out, for a message that
% must stay readable when the list is long.

v = string(v);
if numel(v) <= n
    s = v;
    return;
end
s = [v(1:n); sprintf('... and %d more', numel(v) - n)];
end
