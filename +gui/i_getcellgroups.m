function [labels, groupBy] = i_getcellgroups(sce, parentfig, opts)
%I_GETCELLGROUPS Cell group labels for an analysis, with the usual fallbacks.
%
%   [labels, groupBy] = gui.i_getcellgroups(sce, parentfig)
%   [labels, groupBy] = gui.i_getcellgroups(sce, parentfig, MinGroups=3)
%   [labels, groupBy] = gui.i_getcellgroups(sce, parentfig, Ask=true)
%
%   Resolves the cell grouping a group-wise analysis needs, the way the
%   cell-cell communication callbacks already do it by hand: use
%   SCE.C_CELL_TYPE_TX when it is populated, offer SCE.C_CLUSTER_ID when it
%   is not, and warn when neither can supply enough groups. The dialogs are
%   raised here, so a caller only has to check for an empty return.
%
%   Nothing is written back to SCE. The cluster-id fallback returns
%   "Group<id>" labels rather than assigning them to SCE.C_CELL_TYPE_TX, so
%   running an analysis cannot silently annotate the user's data.
%
%   INPUTS:
%     sce        - SingleCellExperiment
%     parentfig  - parent figure for the dialogs ([] for none)
%     MinGroups  - (2) distinct labels required before a source is usable
%     Ask        - (false) always ask which attribute to group by, via
%                  GUI.I_SELECT1CLASS, instead of preferring cell type
%     Prompt     - prompt for that dialog
%     Notify     - (true) raise the "not enough groups" warning dialog.
%                  Pass false to get the decision without the dialog: the
%                  return value is identical either way.
%
%   NOTIFY=FALSE EXISTS FOR TESTS, AND IT IS NOT A CONVENIENCE. The warning
%   goes through GUI.MYWARNDLG, which is modal - GUI.MYDLG calls
%   WAITFOR(MSGBOX(...)) and blocks until a human clicks OK. A test that
%   reaches the warning path therefore stops the whole suite until somebody
%   notices. GUIGETCELLGROUPSTEST used to do exactly that, and four of its
%   eight cases sat on that dialog.
%
%   OUTPUTS:
%     labels   - NumCells-by-1 string of group labels, or a 0-by-1 string
%                when the user cancelled or no source was usable. Cells
%                with no label are "" and do not count towards MinGroups,
%                but are kept so LABELS stays aligned with the columns of
%                SCE.X.
%     groupBy  - "celltype", "cluster", or the attribute name chosen under
%                Ask=true; "" alongside an empty LABELS.
%
% See also GUI.I_SELECT1CLASS, GUI.CALLBACK_CAUSALCCCNETWORKVIEW.

arguments
    sce
    parentfig = []
    opts.MinGroups (1,1) double {mustBeInteger, mustBePositive} = 2
    opts.Ask (1,1) logical = false
    opts.Prompt (1,1) string = "Select grouping variable:"
    opts.Notify (1,1) logical = true
end

% One place to route the warning through, so a future branch cannot
% reintroduce a blocking dialog by calling MYWARNDLG directly.
warnfun = @(msg) i_notify(opts.Notify, parentfig, msg);

labels = strings(0, 1);
groupBy = "";

if opts.Ask
    [thisc, clabel] = gui.i_select1class(sce, opts.MinGroups <= 1, ...
        char(opts.Prompt), '', parentfig);
    if isempty(thisc), return; end            % cancelled - not an error
    cand = in_tolabels(thisc);
    if in_ngroups(cand) < opts.MinGroups
        warnfun(sprintf( ...
            'Need at least %d cell groups, but "%s" defines %d.', ...
            opts.MinGroups, clabel, in_ngroups(cand)));
        return;
    end
    labels = cand;
    groupBy = string(clabel);
    return;
end

cand = in_tolabels(sce.c_cell_type_tx);
if in_ngroups(cand) >= opts.MinGroups
    labels = cand;
    groupBy = "celltype";
    return;
end

cand = in_tolabels(sce.c_cluster_id);
if in_ngroups(cand) >= opts.MinGroups
    answer = gui.myQuestdlg(parentfig, sprintf( ...
        ['Cell type (C_CELL_TYPE_TX) is undefined.\nWould you like to ', ...
         'use cluster id (C_CLUSTER_ID) to define cell groups?']));
    if ~strcmp(answer, 'Yes'), return; end
    labels = "Group" + cand;
    groupBy = "cluster";
    return;
end

warnfun(sprintf( ...
    'Need at least %d cell groups (SCE.C_CELL_TYPE_TX or SCE.C_CLUSTER_ID).', ...
    opts.MinGroups));
end


function i_notify(on, parentfig, msg)
if on
    gui.myWarndlg(parentfig, msg);
end
end


function s = in_tolabels(c)
% Any of the c_* attributes (cellstr, string, categorical, numeric) as a
% column of strings, with missing values flattened to "" so they read as
% unlabelled rather than as a group of their own.

if isempty(c)
    s = strings(0, 1);
    return;
end
s = string(c);
s = reshape(s, [], 1);
s(ismissing(s)) = "";
s = strip(s);
end


function n = in_ngroups(s)
% Distinct labels, not counting the unlabelled cells.

n = numel(unique(s(s ~= "")));
end
