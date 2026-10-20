function groupIdx = i_pickgroupstolabel(cL, c, parentfig, maxShown)
%I_PICKGROUPSTOLABEL  Choose which groups to label when there are too many.
%
%   groupIdx = gui.i_pickgroupstolabel(cL, c, parentfig)
%   groupIdx = gui.i_pickgroupstolabel(cL, c, parentfig, maxShown)
%
% Called by Label Cell Groups when the plot has more groups than can be
% labelled legibly. Asks whether to label a random MAXSHOWN of them
% (default 50) or to pick the ones to label from a list, which shows each
% group with its cell count.
%
% cL         group names, one per group
% c          group index of every cell, 1..numel(cL)
% parentfig  figure to center the dialogs on
%
% groupIdx is a row vector of indices into cL, sorted, or [] when the user
% cancelled or picked nothing.
%
% See also GUI.MYLISTDLG, GUI.MYQUESTDLG.

if nargin < 4 || isempty(maxShown), maxShown = 50; end
groupIdx = [];
numGroups = numel(cL);
maxShown = min(maxShown, numGroups);

randomChoice = sprintf('Random %d', maxShown);
answer = gui.myQuestdlg(parentfig, sprintf(['There are %d groups, too many ' ...
    'to label at once. Label a random %d of them, or pick which to label?'], ...
    numGroups, maxShown), 'Label Cell Groups', ...
    {randomChoice, 'Pick...', 'Cancel'}, randomChoice);

switch answer
    case randomChoice
        groupIdx = sort(randperm(numGroups, maxShown));
    case 'Pick...'
        cellCounts = accumarray(c(:), 1, [numGroups, 1]);
        items = compose("%s (%d cells)", string(cL(:)), cellCounts);
        [indx, tf] = gui.myListdlg(parentfig, items, 'Label Cell Groups', ...
            [], true, true, [], 'Select the groups to label.');
        if tf == 1
            groupIdx = reshape(sort(indx), 1, []);
        end
    otherwise
        % Cancel, or the dialog was dismissed: label nothing.
end
end
