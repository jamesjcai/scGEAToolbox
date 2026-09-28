function obj = clustercells(obj, k, methodid, forced, sx, opts)
%CLUSTERCELLS  Cluster the cells of a SingleCellExperiment.
%
%   sce = sce.clustercells() clusters the Seurat way: Louvain on a
%   shared-nearest-neighbour graph of the cells' principal components
%   (PKG.E_CELLPCS, then SC_LOUVAIN) at resolution 0.8, letting the graph
%   decide the number of clusters. This is method 'louvainpc', the
%   default. sce.clustercells(k) instead tunes the resolution so that K
%   clusters come out, and sce.clustercells(..., Resolution=G) fixes the
%   resolution at G and ignores K. sce.clustercells(..., Genes=G) adds
%   the genes G to the highly variable genes the components are taken
%   from, for markers too sparse to rank as variable.
%
%   After batch correction (GUI.CALLBACK_HARMONY), 'louvainpc' clusters on
%   the corrected components kept in SCE.STRUCT_CELL_REDUCTIONS.HARMONY, as
%   long as they still have one row per cell. Passing GENES asks for a
%   particular gene set, so it recomputes the components from X instead.
%
%   sce.clustercells(k, methodid) picks another method: 'kmeans',
%   'kmedoids', 'spectclust', 'snndpc', 'mbkmeans' and 'louvain' work on
%   the embedding SCE.S; 'sc3', 'simlr', 'soptsc' and 'sinnlrr' work on
%   the expression matrix. 'louvainpc' also works from the expression
%   matrix and needs no embedding. Note 'louvain' is Louvain on the
%   embedding, not on the principal components.
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

arguments
    obj
    k = []
    methodid = []
    forced = []
    sx = []
    opts.Resolution double = []
    opts.Genes = []
end

if isempty(forced), forced = false; end
if isempty(methodid), methodid = 'louvainpc'; end
methodid = char(methodid);

graphMethods = {'louvainpc'};
% LOUVAINPC lets the graph set the cluster count, so it gets no default K.
if isempty(k) && ~ismember(methodid, graphMethods)
    k = round(obj.NumCells/100, -1);
    if k == 0, k = 1; end
end

embeddingMethods = {'kmeans', 'kmedoids', 'spectclust', 'snndpc', ...
    'mbkmeans', 'louvain'};
expressionMethods = {'sc3', 'simlr', 'soptsc', 'sinnlrr'};

if ~ismember(methodid, [graphMethods, embeddingMethods, expressionMethods])
    error('SingleCellExperiment:clustercells:unknownMethod', ...
        ['Unknown clustering method ''%s''. Principal-component graph: %s. ', ...
        'Embedding-based: %s. Expression-based: %s.'], methodid, ...
        strjoin(graphMethods, ', '), strjoin(embeddingMethods, ', '), ...
        strjoin(expressionMethods, ', '));
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

if ismember(methodid, graphMethods)
    pcs = i_correctedpcs(obj);
    if isempty(pcs) || ~isempty(opts.Genes)
        pcs = pkg.e_cellpcs(obj.X, obj.g, Whitelist=opts.Genes);
    end
    id = sc_louvain(pcs, k, Resolution=opts.Resolution);
elseif usesEmbedding
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

function pcs = i_correctedpcs(obj)
% Batch-corrected components, or [] when there are none or they no longer
% have one row per cell. SELECTCELLS and REMOVECELLS subset them with the
% cells; the row check is for anything that changes the cells without
% going through those.
pcs = [];
r = obj.struct_cell_reductions;
if isstruct(r) && isfield(r, 'harmony') && ~isempty(r.harmony) && ...
        size(r.harmony, 1) == obj.NumCells
    pcs = r.harmony;
end
end
