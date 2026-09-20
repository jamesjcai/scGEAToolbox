function [picked, levels] = i_selectgroupsubset(thisc, clabel, parentfig, sortby)
%I_SELECTGROUPSUBSET Ask which levels of a grouping variable to keep.
%
%   [picked, levels] = gui.i_selectgroupsubset(thisc, clabel, parentfig)
%   [picked, levels] = gui.i_selectgroupsubset(thisc, clabel, parentfig, sortby)
%
%   Inputs:
%     thisc     - per-cell group label vector (string, cellstr, or numeric).
%                 Composite labels from gui.i_selectnclass ("Macrophages |
%                 IL") need no special handling; they are just labels here.
%     clabel    - name of the grouping variable, used in the prompt.
%                 Optional; pass '' or omit to use a generic prompt.
%     parentfig - figure the dialog should be parented to. Optional.
%     sortby    - order the levels are listed in. Optional, one of:
%                   "natural" (default) - natural-sort order, i.e. alphabetic
%                                         with embedded numbers read as numbers
%                   "count"             - most cells first
%                   "none"              - order the levels first appear in
%                   "ask"               - put the choice to the user first,
%                                         which is what a long or lopsided
%                                         annotation needs: the types worth
%                                         keeping are usually the populous
%                                         ones, and alphabetic order buries
%                                         them among the rare ones.
%
%   Outputs:
%     picked - logical column, numel(thisc) long, true for cells whose group
%              was selected. Empty when the user cancelled - callers should
%              treat that as "abort", not as "nothing selected".
%     levels - the selected level names, in the order shown to the user.
%              Bare names, without the cell count the list displays.
%
%   Every caller that lets the user narrow a comparison to a subset of
%   groups goes through here, so the dialog behaves the same everywhere.
%   gui.callback_Violinplot used to call gui.i_selmultidialog instead, which
%   always opens the two-pane shuttle chooser - heavier than the job needs
%   for the handful of groups a violin plot usually compares.
%
%   Each level is listed with its cell count, so a level being tiny - the
%   usual reason to leave one out - is visible in the dialog rather than
%   something the user has to already know.

if nargin < 4 || isempty(sortby), sortby = "natural"; end
if nargin < 3, parentfig = []; end
if nargin < 2, clabel = ''; end

sortby = string(sortby);
picked = [];
levels = strings(0, 1);

thisc = string(thisc);
[ci, cLi] = findgroups(thisc(:));

if sortby == "ask"
    sortby = gui.i_askgrouporder(parentfig);
    if sortby == "", return; end
end

[names, counts] = gui.i_ordergrouplevels(thisc, sortby);
listitems = names + " (n=" + string(counts) + ")";
n = numel(listitems);

if strlength(string(clabel)) > 0
    prompt = sprintf('Select groups to compare (%s):', clabel);
else
    prompt = 'Select groups to compare:';
end

% Preselecting every group is the helpful default for a short list and the
% wrong one for a long one: crossing grouping variables multiplies the level
% count, and past a couple of dozen the user's job becomes deselecting
% rather than selecting.
if n <= 20
    prefersel = listitems;
    initval = 1:n;
else
    prefersel = [];
    initval = 1;   % listdlg wants a valid index, not []
end

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, listitems, prompt, prefersel, true);
else
    [indx, tf] = listdlg('PromptString', {prompt}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listitems, ...
        'InitialValue', initval, ...
        'ListSize', [220, 300]);
end

if tf ~= 1 || isempty(indx), return; end

levels = names(indx);
[found, idx1] = ismember(levels, cLi);
assert(all(found), ...
    'i_selectgroupsubset: a selected level is missing from the group list.');
picked = ismember(ci, idx1);
picked = picked(:);
end
