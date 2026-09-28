function [labels, info] = e_mergesubclusters(labels, targets, opts)
%E_MERGESUBCLUSTERS  Merge groups of cells back into one, the inverse of PKG.E_SUBCLUSTER.
%
%   labels = pkg.e_mergesubclusters(labels, targets) gives every cell whose
%   entry in LABELS is one of TARGETS the same label. Every other cell
%   keeps its label.
%
%   The merged label depends on the class of LABELS:
%     numeric - min(TARGETS). The IDs are then renumbered 1..K in their
%               existing order, so the IDs PKG.E_SUBCLUSTER handed out
%               after max(LABELS) do not leave gaps behind. Undoing a
%               subcluster renumbers nothing else; merging two original
%               clusters shifts the IDs above them down.
%     text    - the common parent when every target is PARENT_{k}, the
%               form PKG.E_SUBCLUSTER writes, so "T cells_{1}" and
%               "T cells_{2}" merge to "T cells".
%
%   labels = pkg.e_mergesubclusters(..., NewLabel=L) names the merged
%   group L instead. It is required for text targets with no common
%   PARENT_{k} form. Numeric IDs are still renumbered afterwards.
%
%   [labels, info] = pkg.e_mergesubclusters(...) also returns the number
%   of cells merged, the merged label, and whether any label outside
%   TARGETS changed in the renumbering.
%
%   See also PKG.E_SUBCLUSTER, PKG.I_SPLITSUBTYPELABEL,
%   GUI.CALLBACK_MERGESUBCLUSTERS.

arguments
    labels (:, 1)
    targets (:, 1)
    opts.NewLabel = []
end

isNumeric = isnumeric(labels);
if ~isNumeric
    labels = string(labels);
    targets = string(targets);
end
targets = unique(targets);
if numel(targets) < 2
    error('pkg:e_mergesubclusters:tooFewTargets', ...
        'Pick at least two groups to merge.');
end
missing = ~ismember(targets, labels);
if any(missing)
    error('pkg:e_mergesubclusters:unknownTarget', ...
        'Group "%s" is not in LABELS. Pick groups the cells carry.', ...
        string(targets(find(missing, 1))));
end

inTargets = ismember(labels, targets);
if ~isempty(opts.NewLabel)
    newLabel = opts.NewLabel;
    if ~isNumeric, newLabel = string(newLabel); end
elseif isNumeric
    newLabel = min(targets);
else
    [parent, suffix, form] = pkg.i_splitsubtypelabel(targets);
    isSubcluster = form == "brace" & ~isnan(str2double(suffix));
    if ~all(isSubcluster) || any(parent ~= parent(1))
        error('pkg:e_mergesubclusters:noCommonParent', ...
            ['The groups do not share a PARENT_{k} name to merge back ' ...
            'into. Pass NewLabel to name the merged group.']);
    end
    newLabel = parent(1);
end

before = labels;
labels(inTargets) = newLabel;

outsideChanged = false;
if isNumeric
    [~, ~, renumbered] = unique(labels);
    renumbered = cast(renumbered, 'like', labels);
    outsideChanged = any(renumbered(~inTargets) ~= before(~inTargets));
    newLabel = renumbered(find(inTargets, 1));
    labels = renumbered;
end

info.NumCells = nnz(inTargets);
info.NewLabel = newLabel;
info.OtherLabelsChanged = outsideChanged;
end
