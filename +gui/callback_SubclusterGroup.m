function [requirerefresh, colorby] = callback_SubclusterGroup(src, ~)
%CALLBACK_SUBCLUSTERGROUP Split one cluster or cell type into subclusters.
%
%   The GUI side of PKG.E_SUBCLUSTER, and the counterpart of Seurat's
%   FindSubCluster: pick a cluster or a cell type, and its cells are
%   re-clustered on their own while every other cell keeps its label.
%   Before this the only ways to split one group were re-clustering every
%   cell, which relabels the groups that were fine, or brushing by hand.
%
%   COLORBY is "cluster" or "celltype", whichever labels were split, so the
%   caller can colour the plot by them.
%
%   See also PKG.E_SUBCLUSTER, GUI.CALLBACK_RECLUSTERCELLS.

requirerefresh = false;
colorby = "cluster";
defaultResolution = 0.5;

[FigureHandle, sce] = gui.gui_getfigsce(src);

hasClusters = numel(unique(sce.c_cluster_id)) > 1;
hasCellTypes = ~isempty(sce.c_cell_type_tx) && ...
    ~all(ismember(string(sce.c_cell_type_tx), ["undetermined", ""]));

if hasClusters && hasCellTypes
    answer = gui.myQuestdlg(FigureHandle, ...
        'Subcluster one cluster or one cell type?', '', ...
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
        'none has a cell type, so there is no group to subcluster. ' ...
        'Cluster the cells first.'], '');
    return;
end

if colorby == "celltype"
    labels = string(sce.c_cell_type_tx(:));
else
    labels = sce.c_cluster_id(:);
end

[grp, grpL] = findgroups(labels);
counts = accumarray(grp, 1);
items = string(grpL) + " (" + counts + " cells)";
[indx, tf] = gui.myListdlg(FigureHandle, items, 'Subcluster Group', ...
    [], false, true, [300, 450], 'Select the group to split into subclusters.');
if tf ~= 1 || isempty(indx), return; end
target = grpL(indx);

resolution = gui.i_askresolution(FigureHandle, defaultResolution, ...
    'Louvain resolution (> 0; Seurat FindSubCluster default 0.5):');
if isempty(resolution), return; end

fw = gui.myWaitbar(FigureHandle);
try
    [labels, info] = pkg.e_subcluster(sce, labels, target, ...
        Resolution=resolution);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if info.NumSubclusters < 2
    gui.myHelpdlg(FigureHandle, sprintf(['Found no subclusters in ' ...
        '"%s" at resolution %g. Try a larger resolution.'], ...
        string(target), resolution));
    return;
end

if colorby == "celltype"
    sce.c_cell_type_tx = labels;
else
    sce.c_cluster_id = labels;
end
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;

gui.myHelpdlg(FigureHandle, sprintf( ...
    '"%s" (%d cells) was split into %d subclusters: %s.', ...
    string(target), info.NumCells, info.NumSubclusters, ...
    strjoin(string(info.NewLabels), ', ')));
end
