function [labels, info] = e_subcluster(sce, labels, target, opts)
%E_SUBCLUSTER  Split one group of cells into subclusters, as Seurat's FindSubCluster.
%
%   labels = pkg.e_subcluster(sce, labels, target) re-clusters the cells of
%   SCE whose entry in LABELS equals TARGET, and returns LABELS with those
%   cells relabelled. Every other cell keeps its label.
%
%   The cells are clustered on their own with SCE.CLUSTERCELLS method
%   'louvainpc': Louvain on a shared-nearest-neighbour graph of principal
%   components. Batch-corrected components in SCE.STRUCT_CELL_REDUCTIONS
%   are used when present, subset to these cells; otherwise the components
%   are recomputed from these cells' counts, so the HVGs are the ones that
%   vary within the group.
%
%   How the new labels are written depends on the class of LABELS:
%     numeric - subcluster 1 keeps TARGET, the others get new IDs after
%               max(LABELS), so every ID stays a number and unique.
%     text    - TARGET_{1}, TARGET_{2}, ..., the suffix cell type
%               annotation already uses for clusters of the same type.
%
%   labels = pkg.e_subcluster(..., Resolution=G) sets the Louvain
%   resolution (default 0.5, FindSubCluster's). MinCells=N (default 20)
%   is the smallest group this will split.
%
%   [labels, info] = pkg.e_subcluster(...) also returns the number of
%   cells split, the subcluster count and the labels given out.
%
%   See also SINGLECELLEXPERIMENT/CLUSTERCELLS, GUI.CALLBACK_SUBCLUSTERGROUP.

arguments
    sce SingleCellExperiment
    labels (:, 1)
    target
    opts.Resolution (1, 1) double {mustBePositive} = 0.5
    opts.MinCells (1, 1) double {mustBePositive, mustBeInteger} = 20
end

if numel(labels) ~= sce.NumCells
    error('pkg:e_subcluster:sizeMismatch', ...
        'LABELS has %d entries but SCE has %d cells. Pass one label per cell.', ...
        numel(labels), sce.NumCells);
end

isNumeric = isnumeric(labels);
if isNumeric
    inGroup = labels == target;
else
    labels = string(labels);
    target = string(target);
    inGroup = labels == target;
end
numInGroup = nnz(inGroup);
if numInGroup < opts.MinCells
    error('pkg:e_subcluster:tooFewCells', ...
        ['Group "%s" has %d cells, fewer than the %d needed to find ', ...
        'subclusters. Pick a larger group.'], string(target), numInGroup, ...
        opts.MinCells);
end

sub = copy(sce);
sub = sub.selectcells(find(inGroup));
sub = sub.clustercells([], 'louvainpc', true, [], Resolution=opts.Resolution);
subId = sub.c_cluster_id(:);
numSub = max(subId);

if isNumeric
    newIds = [target; max(labels) + (1:numSub-1)'];
    newIds = cast(newIds, 'like', labels);
else
    newIds = target + "_{" + (1:numSub)' + "}";
end
if numSub > 1
    labels(inGroup) = newIds(subId);
end

info.NumCells = numInGroup;
info.NumSubclusters = numSub;
info.NewLabels = newIds(1:numSub);
end
