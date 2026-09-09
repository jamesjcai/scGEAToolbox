function [tf, names] = e_hasembedding(sce, ndim)
%E_HASEMBEDDING  Has this SingleCellExperiment actually been embedded?
%
%   tf = PKG.E_HASEMBEDDING(sce) is true when SCE carries at least one real
%   cell embedding. [tf, names] = ... also returns the names of the valid
%   ones, as fields of SCE.STRUCT_CELL_EMBEDDINGS.
%
%   PKG.E_HASEMBEDDING(sce, ndim) restricts to embeddings of that
%   dimensionality, matching the '2d'/'3d' tag in the field name.
%
%   WHY THIS IS NOT isempty(sce.s). The constructor fills S with random
%   coordinates when the caller does not supply any --
%   SingleCellExperiment.m:36 is `if nargin < 3 || isempty(s), s =
%   randn(size(X, 2), 3); end` -- so on an object that has never been
%   embedded, isempty(sce.s) is false and all(sce.s(:) == 0) is false too.
%   Both of the obvious tests pass on randn(nCells, 3).
%
%   That mattered in two places. CLI.CMD_CLUSTER guarded its embedding
%   methods with `isempty(sce.s) || all(sce.s(:) == 0)` and raised "No
%   embedding found in the SCE. Run scgea embed first" -- an error that
%   could never fire, so `scgea cluster` clustered random Gaussian
%   coordinates and reported the result as a cell clustering.
%   SingleCellExperiment.clustercells guarded with `isempty(obj.s)` and had
%   the same hole.
%
%   STRUCT_CELL_EMBEDDINGS is the honest record: it is written only by an
%   actual embedding, so a field that is non-empty and has one row per cell
%   means the work was really done. GUI.I_CHECKEXISTINGEMBED has always
%   tested it this way; this function is that logic in a package both
%   +cli and @SingleCellExperiment can reach.
%
%   See also GUI.I_CHECKEXISTINGEMBED, SC_CLUSTER_S,
%   SINGLECELLEXPERIMENT/EMBEDCELLS.

arguments
    sce SingleCellExperiment
    ndim = []
end

names = {};
if isempty(sce.struct_cell_embeddings)
    tf = false;
    return;
end

slist = fieldnames(sce.struct_cell_embeddings);
valid = false(numel(slist), 1);
for k = 1:numel(slist)
    sx = sce.struct_cell_embeddings.(slist{k});
    if ~isempty(sx) && size(sx, 1) == sce.NumCells
        valid(k) = true;
    end
end
names = slist(valid);

if ~isempty(ndim) && ~isempty(names)
    names = names(contains(string(names), sprintf('%dd', ndim)));
end

tf = ~isempty(names);

end
