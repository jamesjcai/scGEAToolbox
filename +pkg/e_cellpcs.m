function score = e_cellpcs(X, g, npc, numhvg, opts)
%E_CELLPCS  Principal-component scores of cells, the space to cluster in.
%
%   score = pkg.e_cellpcs(X, g) returns a cells x 30 matrix of principal
%   component scores computed from the raw count matrix X (genes x cells)
%   and its gene names G.
%
%   score = pkg.e_cellpcs(X, g, npc, numhvg) sets the number of components
%   (default 30) and of highly variable genes they are taken from (default
%   2000).
%
%   score = pkg.e_cellpcs(..., Whitelist=GENES) adds GENES to the HVG set
%   before the PCA, as EMBEDCELLS does with its whitelist. Canonical markers
%   of rare types are often too sparse to rank as highly variable, and this
%   is how a caller that forced them into an embedding keeps them in the
%   clustering too. Names not in G are ignored.
%
%   The steps follow Seurat's NormalizeData, ScaleData and RunPCA: keep the
%   genes PKG.I_SELECTHVGS picks, the same set EMBEDCELLS embeds; normalise
%   by library size and log1p; scale each gene to zero mean and unit
%   variance, clipped at 10; then PCA. Build a neighbour graph on these
%   scores (SC_LOUVAIN) rather than on a t-SNE or UMAP embedding, which
%   distorts cluster sizes and between-cluster distances.
%
%   See also SC_LOUVAIN, PKG.I_SELECTHVGS, GUI.CALLBACK_SINGLECLICKSOLUTION.

arguments
    X {mustBeNumeric}
    g
    npc (1, 1) double {mustBePositive, mustBeInteger} = 30
    numhvg (1, 1) double {mustBePositive, mustBeInteger} = 2000
    opts.Whitelist = []
end

scaleClip = 10;

g = string(g(:));
keep = true(size(X, 1), 1);
if size(X, 1) > numhvg
    keep = pkg.i_selecthvgs(X, g, numhvg);
end
if ~isempty(opts.Whitelist)
    keep = keep | ismember(g, string(opts.Whitelist(:)));
end
X = X(keep, :);

X = log1p(sc_norm(X));
X(isnan(X)) = 0;
Z = full(X).';

sd = std(Z);
sd(sd == 0) = 1;
Z = (Z - mean(Z)) ./ sd;
Z = min(Z, scaleClip);

npc = min([npc, size(Z, 1) - 1, size(Z, 2)]);
[~, score] = pca(Z, NumComponents=npc);
end
