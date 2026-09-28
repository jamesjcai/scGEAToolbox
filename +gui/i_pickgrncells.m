function cellIdx = i_pickgrncells(sce, parentfig, prompt)
%I_PICKGRNCELLS Choose the cell population a gene network is built from.
%   CELLIDX = gui.i_pickgrncells(SCE, parentfig, PROMPT) asks whether to use
%   all cells or the cells in chosen groups of one grouping, and returns a
%   logical index with one element per cell, or [] if the user cancels.
%   PROMPT is the question shown with the Select Cells... / All Cells
%   buttons. When the cells have no grouping, every cell is used without
%   asking.
%
%   See also gui.i_cellgroupings, gui.callback_BuildGeneNetwork,
%   gui.callback_CompareGeneNetwork.

cellIdx = [];
[names, values] = gui.i_cellgroupings(sce);
if isempty(names)
    cellIdx = true(sce.NumCells, 1);
    return;
end

answer = gui.myQuestdlg(parentfig, prompt, 'Cells', ...
    {'Select Cells...', 'All Cells'}, 'Select Cells...');
switch answer
    case 'Select Cells...'
        % handled below
    case 'All Cells'
        cellIdx = true(sce.NumCells, 1);
        return;
    otherwise
        return;
end

if isscalar(names)
    k = 1;
else
    [k, ok] = gui.myListdlg(parentfig, names, 'Group Cells By', 1, ...
        false, true, [300, 300], 'Select the grouping to choose cells from.');
    if ~ok || isempty(k), return; end
end

[groups, labels] = gui.i_grouplabels(values{k});
[sel, ok] = gui.myListdlg(parentfig, labels, 'Select Cells', [], ...
    true, true, [340, 400], sprintf(['Use the cells in which %s ', ...
    'group(s)? Select one or more.'], names(k)));
if ~ok || isempty(sel), return; end
cellIdx = ismember(values{k}, groups(sel));
end
