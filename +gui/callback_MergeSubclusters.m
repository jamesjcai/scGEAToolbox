function [requirerefresh, colorby] = callback_MergeSubclusters(src, ~)
%CALLBACK_MERGESUBCLUSTERS Merge subclusters back into one group.
%
%   The inverse of GUI.CALLBACK_SUBCLUSTERGROUP, and the GUI side of
%   PKG.E_MERGESUBCLUSTERS. Pick clusters or cell types, and their cells
%   are given one label while every other cell keeps its own.
%
%   Cell types: labels of the form PARENT_{k}, which Subcluster Selected
%   Group and cluster annotation write, are offered back as families, so
%   "T cells_{1}" and "T cells_{2}" go back to "T cells".
%
%   Clusters: a numeric ID keeps no record of the cluster it was split
%   from, so the clusters to merge are picked by hand. They take the
%   smallest of their IDs, and the IDs are renumbered 1..K.
%
%   COLORBY is "cluster" or "celltype", whichever labels were merged, so
%   the caller can colour the plot by them.
%
%   See also PKG.E_MERGESUBCLUSTERS, GUI.CALLBACK_SUBCLUSTERGROUP,
%   GUI.CALLBACK_COLLAPSECELLSUBTYPES.

requirerefresh = false;
colorby = "cluster";

[FigureHandle, sce] = gui.gui_getfigsce(src);

hasClusters = numel(unique(sce.c_cluster_id)) > 1;
hasCellTypes = ~isempty(sce.c_cell_type_tx) && ...
    numel(unique(string(sce.c_cell_type_tx))) > 1;

if hasClusters && hasCellTypes
    answer = gui.myQuestdlg(FigureHandle, ...
        'Merge clusters or cell types?', '', ...
        {'Cluster', 'Cell Type', 'Cancel'}, 'Cluster');
    switch answer
        case 'Cluster'
            colorby = "cluster";
        case 'Cell Type'
            colorby = "celltype";
        otherwise
            return;
    end
elseif hasCellTypes
    colorby = "celltype";
elseif ~hasClusters
    gui.myErrordlg(FigureHandle, ['All cells are in one cluster and ' ...
        'carry one cell type, so there is nothing to merge.'], '');
    return;
end

if colorby == "celltype"
    labels = string(sce.c_cell_type_tx(:));
    [targetSets, tf] = in_pickcelltypefamilies(FigureHandle, labels);
else
    labels = sce.c_cluster_id(:);
    [targetSets, tf] = in_pickclusters(FigureHandle, labels);
end
if ~tf, return; end

numCells = 0;
merged = strings(0, 1);
othersChanged = false;
try
    for k = 1:numel(targetSets)
        [labels, info] = pkg.e_mergesubclusters(labels, targetSets{k});
        numCells = numCells + info.NumCells;
        merged(end+1, 1) = string(info.NewLabel); %#ok<AGROW>
        othersChanged = othersChanged || info.OtherLabelsChanged;
    end
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

if colorby == "celltype"
    sce.c_cell_type_tx = reshape(labels, size(sce.c_cell_type_tx));
else
    sce.c_cluster_id = reshape(labels, size(sce.c_cluster_id));
end
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;

msg = sprintf('%s merged into %s.', pkg.i_plural(numCells, 'cell'), ...
    strjoin('"' + merged + '"', ', '));
if othersChanged
    msg = msg + " The other cluster IDs were renumbered to close the gaps.";
end
gui.myHelpdlg(FigureHandle, msg);
end

function [targetSets, tf] = in_pickcelltypefamilies(FigureHandle, labels)
% Families of PARENT_{k} labels, with PARENT itself when cells still carry it.
targetSets = {};
[ulabels, ~, back] = unique(labels);
counts = accumarray(back, 1);
[parent, suffix, form] = pkg.i_splitsubtypelabel(ulabels);
isSubcluster = form == "brace" & ~isnan(str2double(suffix));
parents = unique(parent(isSubcluster));

families = {};
items = strings(0, 1);
for k = 1:numel(parents)
    members = ulabels(isSubcluster & parent == parents(k));
    if ismember(parents(k), ulabels)
        members = [parents(k); members]; %#ok<AGROW>
    end
    if numel(members) < 2, continue; end
    families{end+1} = members; %#ok<AGROW>
    items(end+1, 1) = sprintf('%s  ->  %s   (%s)', ...
        strjoin(members, ', '), parents(k), ...
        pkg.i_plural(sum(counts(ismember(ulabels, members))), 'cell')); %#ok<AGROW>
end

if isempty(families)
    gui.myHelpdlg(FigureHandle, ['No two cell types share a name of the ' ...
        'form "Type_{1}", "Type_{2}", so there are no subclusters to ' ...
        'merge back. Annotate > Collapse Cell Subtypes handles other ' ...
        'subtype labels.']);
    tf = false;
    return;
end

[pick, tf] = gui.myListdlg(FigureHandle, items, 'Merge Subclusters', ...
    1:numel(items), true, true, [460, 320], ...
    'Select the subclusters to merge back into their cell type.');
tf = tf == 1 && ~isempty(pick);
if tf, targetSets = families(pick); end
end

function [targetSets, tf] = in_pickclusters(FigureHandle, labels)
targetSets = {};
[grp, grpL] = findgroups(labels);
counts = accumarray(grp, 1);
items = string(grpL) + " (" + counts + " cells)";
[pick, tf] = gui.myListdlg(FigureHandle, items, 'Merge Subclusters', ...
    [], true, true, [300, 450], ...
    'Select two or more clusters to merge into one.');
tf = tf == 1 && ~isempty(pick);
if ~tf, return; end
if numel(pick) < 2
    gui.myHelpdlg(FigureHandle, 'Select at least two clusters to merge.');
    tf = false;
    return;
end
targetSets = {grpL(pick)};
end
