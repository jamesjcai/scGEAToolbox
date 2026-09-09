function obj = clustercells(obj, k, methodid, forced, sx)
%CLUSTERCELLS  Cluster the cells of a SingleCellExperiment.
%
%   sce = sce.clustercells(k) clusters into K groups using the cell
%   embedding. sce.clustercells(k, methodid) picks the method:
%   'kmeans' (default), 'kmedoids', 'spectclust', 'snndpc' or 'mbkmeans'
%   work on the embedding; 'sc3', 'simlr', 'soptsc' and 'sinnlrr' work on
%   the expression matrix and need no embedding.
%
%   sce.clustercells(k, methodid, forced) recomputes even when a
%   clustering is already present. SX supplies an embedding to use in
%   place of SCE.S.
%
%   THREE THINGS THAT WERE WRONG HERE.
%
%   (1) The guard was `isempty(obj.s)`, which is never true: the
%   constructor fills S with randn(nCells, 3). So an object that had never
%   been embedded passed the check and got clustered on random
%   coordinates. PKG.E_HASEMBEDDING tests STRUCT_CELL_EMBEDDINGS, which
%   only a real embedding writes.
%
%   (2) The recompute guard was `isempty(obj.c_cluster_id) || forced`, and
%   the constructor presets C_CLUSTER_ID to ones(nCells, 1), so that was
%   never empty either. Without FORCED this method silently did nothing:
%   `sce.clustercells(12)` -- the form shown in the help of
%   LLM.E_CELLTYPEANNO and SC_ANNOTATECELLS -- returned one cluster having
%   been asked for twelve, and printed nothing. Every in-tree caller had
%   worked around it by always passing forced = true.
%
%   (3) The switch had no OTHERWISE and covered only four of the nine
%   methods the two clustering functions accept, so 'kmedoids',
%   'spectclust' and 'mbkmeans' fell through with ID unassigned and failed
%   on the next line with "Unable to index into 'id'".
%
%   See also SC_CLUSTER_S, SC_CLUSTER_X, PKG.E_HASEMBEDDING.

if nargin < 5, sx = []; end
if nargin < 4 || isempty(forced), forced = false; end
if nargin < 3 || isempty(methodid), methodid = 'kmeans'; end
if nargin < 2 || isempty(k)
    k = round(obj.NumCells/100, -1);
    if k == 0, k = 1; end
end

embeddingMethods = {'kmeans', 'kmedoids', 'spectclust', 'snndpc', 'mbkmeans'};
expressionMethods = {'sc3', 'simlr', 'soptsc', 'sinnlrr'};

if ~ismember(methodid, [embeddingMethods, expressionMethods])
    error('SingleCellExperiment:clustercells:unknownMethod', ...
        ['Unknown clustering method ''%s''. Embedding-based: %s. ', ...
        'Expression-based: %s.'], methodid, ...
        strjoin(embeddingMethods, ', '), strjoin(expressionMethods, ', '));
end

usesEmbedding = ismember(methodid, embeddingMethods);
if usesEmbedding && isempty(sx) && ~pkg.e_hasembedding(obj)
    error('SingleCellExperiment:clustercells:noEmbedding', ...
        ['Method ''%s'' needs a cell embedding and this object has none. ', ...
        'Run sce = sce.embedcells(...) first, pass one as SX, or use an ', ...
        'expression-based method (%s).'], methodid, ...
        strjoin(expressionMethods, ', '));
end

alreadyClustered = isfield(obj.struct_cell_clusterings, methodid) && ...
    ~isempty(obj.struct_cell_clusterings.(methodid));
if alreadyClustered && ~forced
    return;
end

if usesEmbedding
    if isempty(sx)
        id = sc_cluster_s(obj.s, k, 'type', methodid);
    else
        id = sc_cluster_s(sx, k, 'type', methodid);
    end
else
    id = sc_cluster_x(obj.X, k, 'type', methodid);
end

obj.c_cluster_id = id(:);
obj.struct_cell_clusterings.(methodid) = obj.c_cluster_id;
end
