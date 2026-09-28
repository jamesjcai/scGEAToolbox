function [insce, names] = i_pickworkspacesces(parentfig, prompt, exclude)
% I_PICKWORKSPACESCES - Pick one or more SCE variables from the base workspace.
%
%   [insce, names] = gui.i_pickworkspacesces(parentfig, prompt)
%   [insce, names] = gui.i_pickworkspacesces(parentfig, prompt, exclude)
%
% Lists every scalar SingleCellExperiment in the base workspace, in natural
% order, each with its cell and gene counts. INSCE is a cell array of the
% selected objects and NAMES their variable names, as a string array; both
% are empty when nothing qualifies or the dialog is cancelled. How many must
% be selected is the caller's to check.
%
% EXCLUDE is an SCE left out of the list, for a caller that already has it:
% a workspace variable that is the same handle as the dataset open in the
% window would otherwise be offered for merging with itself.
%
% see also: gui.i_mergesces, gui.sc_openscedlg

if nargin < 3, exclude = []; end
insce = {};
names = strings(0);

vars = evalin('base', 'whos');
vars = vars(strcmp({vars.class}, 'SingleCellExperiment'));
vars = vars(arrayfun(@(v) isequal(v.size, [1, 1]), vars));

cand = cell(1, numel(vars));
keep = true(1, numel(vars));
for k = 1:numel(vars)
    cand{k} = evalin('base', vars(k).name);
    keep(k) = isempty(exclude) || cand{k} ~= exclude;
end
cand = cand(keep);
vars = vars(keep);

if isempty(vars)
    msg = 'There are no SCE variables in the base workspace.';
    if ~isempty(exclude)
        msg = ['There are no SCE variables in the base workspace other ' ...
            'than the dataset open in this window.'];
    end
    gui.myWarndlg(parentfig, msg);
    return;
end

[candnames, idx] = natsort(string({vars.name}));
cand = cand(idx);
items = strings(1, numel(cand));
for k = 1:numel(cand)
    items(k) = sprintf("%s (%d cells, %d genes)", candnames(k), ...
        cand{k}.NumCells, cand{k}.NumGenes);
end

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, items, 'SCE Variables', [], ...
        true, true, [480, 300], prompt);
else
    [indx, tf] = listdlg('PromptString', {prompt}, 'ListString', items, ...
        'SelectionMode', 'multiple', 'ListSize', [460, 300]);
end
if tf ~= 1 || isempty(indx), return; end
insce = cand(indx);
names = candnames(indx);
end
