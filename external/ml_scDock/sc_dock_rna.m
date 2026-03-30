function result = sc_dock_rna(X, g, varargin)
% SC_DOCK_RNA  Module 1 of scDock: scRNA-seq preprocessing and annotation.
%
% Performs QC filtering, normalization, HVG selection, PCA (with geometric
% elbow PC selection), optional batch correction, UMAP embedding,
% clustering, and cell-type annotation.
%
% USAGE:
%   result = sc_dock_rna(X, g)
%   result = sc_dock_rna(X, g, 'batch_id', batch_vec, 'species', 'mouse')
%
% INPUTS:
%   X         - genes x cells raw count matrix (sparse or dense)
%   g         - gene list (string array or cell array of chars)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'batch_id'        - 1 x ncells batch label vector ([] = no correction)
%   'species'         - 'human' (default) or 'mouse'
%   'libsz_cutoff'    - minimum library size per cell (default 500)
%   'mt_ratio'        - max mitochondrial read ratio (default 0.15)
%   'min_cells'       - minimum cells a gene must appear in (default 10)
%   'n_hvg'           - number of HVGs to use for PCA (default 2000)
%   'n_pcs'           - max PCs to compute; 0 = auto via elbow (default 0)
%   'n_clusters'      - number of clusters for sc_cluster_s (default 10)
%   'clust_method'    - clustering method for sc_cluster_s (default 'kmeans')
%   'batch_method'    - batch correction method: 'harmony' | 'none'
%                       (default 'harmony' when batch_id provided)
%
% OUTPUT:
%   result  - struct with fields:
%     .X_norm       - normalized log1p expression (genes x cells, HVGs only)
%     .g_hvg        - gene names corresponding to rows of X_norm
%     .pca_scores   - PCA scores (cells x n_pcs)
%     .pca_var      - variance explained per PC
%     .n_pcs        - number of PCs selected by geometric elbow
%     .s_umap       - UMAP embedding (cells x 3)
%     .c_cluster_id - cluster assignments (cells x 1)
%     .c_cell_type  - cell-type annotation strings (cells x 1)
%     .cell_kept    - logical index of cells kept after QC (original order)
%
% DEPENDENCIES (scGEAToolbox_dev):
%   sc_qcfilter, sc_norm, sc_hvg, pkg.e_randPCA,
%   run.ml_Harmony2, sc_umap, sc_cluster_s, sc_celltypeanno
%
% See also: SC_DOCK_CCC, SC_DOCK_VINA, SC_DOCK

p = inputParser;
addRequired(p, 'X', @(x) isnumeric(x) || islogical(x));
addRequired(p, 'g');
addParameter(p, 'batch_id',     [],         @(x) isempty(x) || isvector(x));
addParameter(p, 'species',      'human',    @(x) ismember(x, {'human','mouse'}));
addParameter(p, 'libsz_cutoff', 500,        @(x) isscalar(x) && x >= 0);
addParameter(p, 'mt_ratio',     0.15,       @(x) isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'min_cells',    10,         @(x) isscalar(x) && x >= 1);
addParameter(p, 'n_hvg',        2000,       @(x) isscalar(x) && x >= 1);
addParameter(p, 'n_pcs',        0,          @(x) isscalar(x) && x >= 0);
addParameter(p, 'n_clusters',   10,         @(x) isscalar(x) && x >= 2);
addParameter(p, 'clust_method', 'kmeans',   @ischar);
addParameter(p, 'batch_method', 'harmony',  @(x) ismember(x, {'harmony','none'}));
parse(p, X, g, varargin{:});
opt = p.Results;

if iscellstr(g), g = string(g); end %#ok<ISCLSTR>
n_cells_orig = size(X, 2);

% -------------------------------------------------------------------------
% Step 1: QC filtering
% -------------------------------------------------------------------------
fprintf('[sc_dock_rna] Step 1/7: QC filtering (%d genes x %d cells)...\n', ...
    size(X,1), size(X,2));
[X, g, kept_idx_list] = sc_qcfilter(X, g, opt.libsz_cutoff, ...
    opt.mt_ratio, opt.min_cells);

% Build a single logical index mapping original cells to kept cells
cell_kept = true(n_cells_orig, 1);
for k = 1:numel(kept_idx_list)
    idx = kept_idx_list{k};
    if islogical(idx)
        tmp = cell_kept;
        tmp(cell_kept) = idx;
        cell_kept = tmp;
    else
        tmp = false(sum(cell_kept), 1);
        tmp(idx) = true;
        cell_kept2 = cell_kept;
        cell_kept2(cell_kept) = tmp;
        cell_kept = cell_kept2;
    end
end

if ~isempty(opt.batch_id)
    batch_id = opt.batch_id(cell_kept);
else
    batch_id = [];
end

fprintf('    -> %d genes x %d cells after QC.\n', size(X,1), size(X,2));

% -------------------------------------------------------------------------
% Step 2: Normalization (library-size + log1p)
% -------------------------------------------------------------------------
fprintf('[sc_dock_rna] Step 2/7: Normalizing...\n');
X_norm = sc_norm(X, 'libsize');
X_norm = log1p(X_norm);

% -------------------------------------------------------------------------
% Step 3: HVG selection
% -------------------------------------------------------------------------
fprintf('[sc_dock_rna] Step 3/7: Selecting HVGs...\n');
n_hvg = min(opt.n_hvg, size(X_norm, 1));
[~, X_hvg, g_hvg] = sc_hvg(X_norm, g, true, false, false);
X_hvg = X_hvg(1:n_hvg, :);
g_hvg  = g_hvg(1:n_hvg);
fprintf('    -> %d HVGs selected.\n', n_hvg);

% -------------------------------------------------------------------------
% Step 4: PCA
% -------------------------------------------------------------------------
fprintf('[sc_dock_rna] Step 4/7: Running PCA...\n');
n_pcs_max = min(opt.n_pcs > 0, 1) * opt.n_pcs + ...
            (opt.n_pcs == 0) * min(50, size(X_hvg, 2) - 1);
n_pcs_max = max(n_pcs_max, 10);

% Center rows (genes), then transpose so pca is on cells
Xc = double(X_hvg) - mean(double(X_hvg), 2);  % center genes
[~, S, V] = pkg.e_randPCA(Xc, n_pcs_max);    % V is cells x n_pcs_max
pca_scores  = V;                                % cells x n_pcs_max
pca_var     = diag(S).^2;
pca_var_pct = pca_var ./ sum(pca_var) * 100;

% Select number of PCs using geometric elbow method
if opt.n_pcs == 0
    n_pcs = i_geom_elbow(pca_var_pct);
    fprintf('    -> Geometric elbow selected %d PCs.\n', n_pcs);
else
    n_pcs = opt.n_pcs;
    fprintf('    -> Using %d PCs (user-specified).\n', n_pcs);
end
pca_scores = pca_scores(:, 1:n_pcs);

% -------------------------------------------------------------------------
% Step 5: Batch correction (optional)
% -------------------------------------------------------------------------
if ~isempty(batch_id) && strcmp(opt.batch_method, 'harmony')
    fprintf('[sc_dock_rna] Step 5/7: Harmony batch correction...\n');
    pca_scores = run.ml_Harmony2(pca_scores, batch_id);
else
    fprintf('[sc_dock_rna] Step 5/7: Skipping batch correction.\n');
end

% -------------------------------------------------------------------------
% Step 6: UMAP embedding (on PCA scores)
% -------------------------------------------------------------------------
fprintf('[sc_dock_rna] Step 6/7: UMAP embedding...\n');
% sc_umap normalizes internally; we pass already-normalized PCA scores
% directly by disabling its internal norm/log
s_umap = sc_umap(pca_scores', 3, false, false);

% -------------------------------------------------------------------------
% Step 7: Clustering + cell-type annotation
% -------------------------------------------------------------------------
fprintf('[sc_dock_rna] Step 7/7: Clustering and cell-type annotation...\n');
c_cluster_id = sc_cluster_s(pca_scores, opt.n_clusters, opt.clust_method);
[c_cell_type] = sc_celltypeanno(X_norm, g, c_cluster_id, opt.species);

% -------------------------------------------------------------------------
% Pack results
% -------------------------------------------------------------------------
result.X_norm       = X_hvg;
result.g_hvg        = g_hvg;
result.pca_scores   = pca_scores;
result.pca_var      = pca_var_pct;
result.n_pcs        = n_pcs;
result.s_umap       = s_umap;
result.c_cluster_id = c_cluster_id;
result.c_cell_type  = c_cell_type;
result.cell_kept    = cell_kept;

fprintf('[sc_dock_rna] Done. %d cells, %d clusters.\n', ...
    size(pca_scores,1), numel(unique(c_cluster_id)));
end

% =========================================================================
% Local helper: geometric elbow method for PC selection
% Identifies PC index with maximum perpendicular distance from the line
% connecting the first and last points of the variance-explained curve.
% Ref: Zhuang et al. 2022
% =========================================================================
function n = i_geom_elbow(var_pct)
n_max = numel(var_pct);
% Line endpoints
p1 = [1,          var_pct(1)];
p2 = [n_max,      var_pct(end)];
d  = p2 - p1;
d_norm = d / norm(d);

% Perpendicular distance from each point to the line
dists = zeros(n_max, 1);
for i = 1:n_max
    v = [i, var_pct(i)] - p1;
    dists(i) = abs(v(1)*d_norm(2) - v(2)*d_norm(1));
end
[~, n] = max(dists);
n = max(n, 2);  % always use at least 2 PCs
end
