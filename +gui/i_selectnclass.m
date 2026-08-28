function [thisc, clabel, partlabels] = i_selectnclass(sce, allowsingle, ...
    promptstr, prefersel, parentfig, sep)
%I_SELECTNCLASS  Pick one or more grouping variables and cross them.
%
%   [thisc, clabel, partlabels] = gui.i_selectnclass(sce)
%   [...] = gui.i_selectnclass(sce, allowsingle, promptstr, prefersel, parentfig, sep)
%
% The multi-select sibling of GUI.I_SELECT1CLASS. Where that one returns a
% single attribute as stored, this one crosses every selected variable into
% one composite label per cell - a cell that is a Macrophage from batch IL
% becomes "Macrophages | IL" - so a caller that groups by THISC gets one
% group per observed combination.
%
% Only combinations that actually occur are produced. Crossing is done per
% cell, never as a full product of levels, so empty combinations never
% appear as empty groups.
%
% Outputs
%   thisc       NumCells x 1 string, the composite label of each cell.
%               A single selection gives string(attribute), so the
%               one-variable case behaves as I_SELECT1CLASS does.
%   clabel      the composite variable name, e.g. 'Cell Type x Batch ID',
%               suitable for an axis label or a dialog title.
%   partlabels  1 x k cellstr of the chosen variable names, in list order.
%
% Inputs match I_SELECT1CLASS, plus SEP (default " | "), the separator used
% to join the parts. It is a parameter because a level name that itself
% contains the separator would otherwise be ambiguous.
%
% Returns empty THISC when the user cancels, so the usual caller guard
% `if isempty(thisc), return; end` works unchanged.
%
% See also GUI.I_SELECT1CLASS, GUI.I_SELECTNSTATES.

if nargin < 6 || isempty(sep), sep = " | "; end
if nargin < 5, parentfig = []; end
if nargin < 4, prefersel = ''; end
if nargin < 3 || isempty(promptstr)
    promptstr = 'Select one or more grouping variables:';
end
if nargin < 2 || isempty(allowsingle), allowsingle = true; end
if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end

% Above this many groups the violin/heatmap views stop being readable, so
% ask rather than render something unusable. Crossing two attributes with
% many levels each (cluster id x cell type, say) reaches the thousands.
MAXLEVELS = 60;
MISSINGTAG = "(n/a)";

thisc = [];
clabel = '';
partlabels = {};

[listitems, wsvars, wsinfo] = i_classlistitems(sce, allowsingle);

y = false;
if ~isempty(prefersel)
    [y, idx] = ismember(prefersel, listitems);
end

if gui.i_isuifig(parentfig)
    % allowmulti = true, and passed explicitly. Relying on myListdlg's
    % default is what let a multi-select listbox appear in front of code
    % that could only honour one pick.
    [indx2, tf2] = gui.myListdlg(parentfig, listitems, ...
        promptstr, prefersel, true);
else
    if y
        [indx2, tf2] = listdlg('PromptString', {promptstr}, ...
            'SelectionMode', 'multiple', 'ListString', listitems, ...
            'ListSize', [220 300], 'InitialValue', idx, 'Name', ' ');
    else
        [indx2, tf2] = listdlg('PromptString', {promptstr}, ...
            'SelectionMode', 'multiple', 'ListString', listitems, ...
            'ListSize', [220 300], 'Name', ' ');
    end
end

if tf2 ~= 1
    if ~isempty(parentfig), figure(parentfig); end
    return;
end

% List order, deduplicated: the same two variables must give the same
% composite however the user happened to click them.
indx2 = sort(unique(indx2(:))).';
k = numel(indx2);
if k == 0, return; end

S = strings(sce.NumCells, k);
partlabels = cell(1, k);
for j = 1:k
    [v, partlabels{j}] = i_getone(indx2(j));
    if isempty(v)                   % cancelled the workspace-variable pick
        thisc = []; clabel = ''; partlabels = {};
        return;
    end
    S(:, j) = string(v(:));
end

thisc = gui.i_crossgroups(S, sep, MISSINGTAG);
clabel = char(strjoin(string(partlabels), ' x '));

nlev = numel(unique(thisc));
if nlev > MAXLEVELS
    answer = gui.myQuestdlg(parentfig, sprintf(['Grouping by %s gives %d ' ...
        'groups. Plots with that many groups are usually unreadable. ' ...
        'Continue?'], clabel, nlev), '', {'Yes', 'No'}, 'No', 'warning');
    if ~strcmp(answer, 'Yes')
        thisc = []; clabel = ''; partlabels = {};
        return;
    end
end


    function [c, thislabel] = i_getone(indx)
        % One selected variable, resolved the same way i_select1class does.
        c = [];
        thislabel = listitems{indx};
        switch thislabel
            case 'Current Class (C)'
                c = sce.c;
            case 'Cluster ID'
                c = sce.c_cluster_id;
            case 'Batch ID'
                c = sce.c_batch_id;
            case 'Cell Type'
                c = sce.c_cell_type_tx;
            case 'Cell Cycle Phase'
                c = sce.c_cell_cycle_tx;
            case 'Workspace Variable...'
                if gui.i_isuifig(parentfig)
                    [wi, wtf] = gui.myListdlg(parentfig, wsinfo(1, :), ...
                        'Select variable:', [], false);
                else
                    [wi, wtf] = listdlg('PromptString', {'Select variable:'}, ...
                        'liststring', wsinfo(1, :), 'SelectionMode', 'single');
                end
                if wtf == 1
                    c = evalin('base', wsvars(wi).name);
                    thislabel = wsvars(wi).name;
                end
        end
    end

end
