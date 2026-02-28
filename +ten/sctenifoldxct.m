function [T] = sctenifoldxct(sce_ori, celltype1, celltype2, twosided, varargin)
%SCTENIFOLDXCT  Native MATLAB implementation of scTenifoldXct.
%   Detects ligand-receptor-mediated cell-cell interactions by spectral
%   manifold alignment of gene regulatory networks.  No Python required.
%
%   T = ten.sctenifoldxct(SCE, CELLTYPE1, CELLTYPE2)
%   T = ten.sctenifoldxct(SCE, CELLTYPE1, CELLTYPE2, TWOSIDED)
%   T = ten.sctenifoldxct(SCE, CELLTYPE1, CELLTYPE2, TWOSIDED, Name, Value)
%
%   Inputs:
%     SCE       - SingleCellExperiment object
%     CELLTYPE1 - source cell type label (string)
%     CELLTYPE2 - target cell type label (string)
%     TWOSIDED  - logical, run both directions (default: true)
%
%   Name-Value pairs:
%     'n_dim'   - spectral embedding dimensions (default: 50)
%     'mu'      - cross-cell-type correspondence weight (default: 0.9)
%     'pval'    - p-value threshold for null test, 1.0 returns all pairs
%                 (default: 1.0)
%     'verbose' - print progress messages (default: true)
%
%   Output:
%     T - table with columns: ligand, receptor, dist, p_value
%         If twosided=true, T is a cell array {T1, T2} where T1 is
%         CELLTYPE1→CELLTYPE2 and T2 is CELLTYPE2→CELLTYPE1.
%
%   Algorithm:
%     Builds partial-correlation GRNs for each cell type via sc_pcnetpar,
%     then performs spectral manifold alignment via graph Laplacian
%     eigenvectors on the combined [GRN_source, W_LR; W_LR', GRN_target]
%     weight matrix.  Ligand-receptor pairs from the built-in database act
%     as correspondences that pull ligand and receptor gene embeddings
%     together.  Significance is assessed by a nonparametric left-tail null
%     test comparing candidate L-R distances to a sampled background.
%
%   This function is a drop-in MATLAB replacement for +run/py_scTenifoldXct.
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742

if nargin < 4, twosided = true; end

p = inputParser;
addOptional(p, 'n_dim',   50,   @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',      0.9,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'pval',    1.0,  @(x) isnumeric(x) && x >= 0);
addOptional(p, 'verbose', true, @islogical);
parse(p, varargin{:});

n_dim   = round(p.Results.n_dim);
mu_     = p.Results.mu;
pval_t  = p.Results.pval;
verbose = p.Results.verbose;

% ── 1. Prepare expression matrices ───────────────────────────────────────
sce = copy(sce_ori);
idx = sce.c_cell_type_tx == celltype1 | sce.c_cell_type_tx == celltype2;
sce = sce.selectcells(idx);

is_s = sce.c_cell_type_tx == celltype1;
is_t = sce.c_cell_type_tx == celltype2;

X = sce.X;
if issparse(X), X = full(X); end
X   = single(X);
g   = sce.g;
X_s = X(:, is_s);   % genes × source cells
X_t = X(:, is_t);   % genes × target cells

if verbose
    fprintf('[sctenifoldxct] %s: %d genes × %d cells\n', ...
        celltype1, size(X_s,1), size(X_s,2));
    fprintf('[sctenifoldxct] %s: %d genes × %d cells\n', ...
        celltype2, size(X_t,1), size(X_t,2));
end

% ── 2. Load ligand-receptor database ─────────────────────────────────────
pw1     = fileparts(mfilename('fullpath'));
lr_mat  = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.mat');
lr_txt  = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.txt');

if exist(lr_mat, 'file')
    db     = load(lr_mat, 'ligand', 'receptor');
    lig_db = upper(string(db.ligand(:)));
    rec_db = upper(string(db.receptor(:)));
else
    % Fallback: read from text file (tab-delimited, cols 3 and 5)
    T_lr   = readtable(lr_txt, 'FileType', 'text', 'Delimiter', '\t');
    lig_db = upper(string(T_lr{:,3}));
    rec_db = upper(string(T_lr{:,5}));
end
if verbose
    fprintf('[sctenifoldxct] L-R database: %d pairs loaded.\n', numel(lig_db));
end

% ── 3. Build GRNs ────────────────────────────────────────────────────────
if verbose, fprintf('[sctenifoldxct] Building GRN: %s ...\n', celltype1); end
A_s = sc_pcnetpar(X_s);
A_s = A_s ./ max(abs(A_s(:)));
A_s = ten.e_filtadjc(A_s, 0.75, false);   % dense, thresholded

if verbose, fprintf('[sctenifoldxct] Building GRN: %s ...\n', celltype2); end
A_t = sc_pcnetpar(X_t);
A_t = A_t ./ max(abs(A_t(:)));
A_t = ten.e_filtadjc(A_t, 0.75, false);

% ── 4. Manifold alignment + null test ────────────────────────────────────
if verbose
    fprintf('[sctenifoldxct] Aligning %s → %s ...\n', celltype1, celltype2);
end
T1 = i_xct(X_s, X_t, g, A_s, A_t, lig_db, rec_db, n_dim, mu_, pval_t, verbose);

if twosided
    if verbose
        fprintf('[sctenifoldxct] Aligning %s → %s ...\n', celltype2, celltype1);
    end
    T2 = i_xct(X_t, X_s, g, A_t, A_s, lig_db, rec_db, n_dim, mu_, pval_t, verbose);
    T  = {T1, T2};
else
    T = T1;
end

end % main function


%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────

function T = i_xct(X_s, X_t, g, A_s, A_t, lig_db, rec_db, n_dim, mu_, pval_t, verbose)
%I_XCT  Spectral manifold alignment for one direction (source → target).

ng   = size(X_s, 1);   % number of genes (same for source and target)
g_up = upper(string(g(:)));

% ── Build sparse L-R correspondence matrix W12 (ng × ng) ─────────────────
n_lr = numel(lig_db);
li_idx = zeros(n_lr, 1);
ri_idx = zeros(n_lr, 1);
n_valid = 0;
for k = 1:n_lr
    li = find(g_up == lig_db(k), 1);
    ri = find(g_up == rec_db(k), 1);
    if ~isempty(li) && ~isempty(ri)
        n_valid = n_valid + 1;
        li_idx(n_valid) = li;
        ri_idx(n_valid) = ri;
    end
end
li_idx = li_idx(1:n_valid);
ri_idx = ri_idx(1:n_valid);

if n_valid == 0
    warning('sctenifoldxct:noPairs', ...
        'No L-R pairs found in gene list. Returning empty table.');
    T = table();
    return;
end
if verbose
    fprintf('[sctenifoldxct]   %d L-R pairs matched in data.\n', n_valid);
end

% Aggregate duplicate L-R pairs (sum weights)
W12 = sparse(li_idx, ri_idx, ones(n_valid,1), ng, ng);

% ── Build symmetric GRN weight matrices ──────────────────────────────────
% Symmetrize and shift by +1 on diagonal for graph connectivity
%   (follows convention in ten.i_ma used by sctenifoldnet)
W11 = sparse(0.5*(A_s + A_s'));
W22 = sparse(0.5*(A_t + A_t'));

% Scale mu to balance within/across-type connection strengths
w11_sum = sum(abs(W11(:)));
w22_sum = sum(abs(W22(:)));
w12_sum = sum(abs(W12(:)));
if w12_sum == 0
    mu_scale = 0;
else
    mu_scale = mu_ * (w11_sum + w22_sum) / (2 * w12_sum);
end

% Block weight matrix: [source | cross; cross' | target]
W = [W11,                mu_scale .* W12;
     mu_scale .* W12',   W22            ];

% ── Graph Laplacian (unnormalized) ────────────────────────────────────────
d = full(sum(abs(W), 2));   % degree vector (using |W| for signed graphs)
L = spdiags(d, 0, 2*ng, 2*ng) - W;

% ── Spectral embedding via n_dim smallest non-trivial eigenvectors ────────
n_ev   = min(n_dim + 4, 2*ng - 2);   % request a few extra for filtering
opts.isreal  = true;
opts.issym   = true;
opts.tol     = 1e-6;
[V, D_ev] = eigs(L, n_ev, 'smallestreal', opts);
ev = real(diag(D_ev));
V  = real(V);
[ev, ord] = sort(ev, 'ascend');
V  = V(:, ord);

% Discard trivial (near-zero, constant) eigenvectors
keep = ev >= 1e-8;
V    = V(:, keep);
ev   = ev(keep);   %#ok<NASGU>

if size(V, 2) < n_dim
    n_dim_used = size(V, 2);
    if verbose
        warning('sctenifoldxct:dimReduced', ...
            'Only %d non-trivial eigenvectors available; n_dim reduced from %d.', ...
            n_dim_used, n_dim);
    end
else
    n_dim_used = n_dim;
end
V = V(:, 1:n_dim_used);   % (2*ng) × n_dim_used

P_s = V(1:ng,    :);   % source gene embeddings  (ng × n_dim)
P_t = V(ng+1:end,:);   % target gene embeddings  (ng × n_dim)

% ── Candidate L-R pair distances ─────────────────────────────────────────
n_cand  = n_valid;
cand_d  = zeros(n_cand, 1, 'single');
for k = 1:n_cand
    diff_k   = P_s(li_idx(k),:) - P_t(ri_idx(k),:);
    cand_d(k) = sqrt(sum(diff_k.^2));
end

% ── Null distribution: sample random non-L-R gene pairs ──────────────────
% L-R pairs are ~n_valid out of ng^2 total; random sampling rarely hits them
n_null = max(50000, 100 * n_cand);
n_null = min(n_null, ng^2 - n_cand);   % cap at available pairs

rng_state = rng;   % save rng state (non-destructive)
rand_i = randi(ng, n_null, 1);
rand_j = randi(ng, n_null, 1);
rng(rng_state);

diff_null = P_s(rand_i,:) - P_t(rand_j,:);
null_d    = sqrt(sum(diff_null.^2, 2));

% ── Left-tail null test ───────────────────────────────────────────────────
%   p_value(k) = fraction of null distances ≤ cand_d(k)
%   Small p_value → candidate distance unusually small → strong interaction
p_vals = sum(null_d(:) <= cand_d(:)', 1)' ./ numel(null_d);

% ── Assemble output table ─────────────────────────────────────────────────
lig_out  = g(li_idx);
rec_out  = g(ri_idx);

% Lookup per-pair W12 weights (number of duplicate LR entries, if any)
w12_vals = full(W12(sub2ind([ng,ng], li_idx, ri_idx)));

T = table(lig_out, rec_out, double(cand_d), w12_vals, p_vals, ...
    'VariableNames', {'ligand','receptor','dist','correspondence','p_value'});
T = T(T.p_value <= pval_t, :);
T = sortrows(T, 'dist', 'ascend');

if verbose
    fprintf('[sctenifoldxct]   %d significant L-R pairs (pval ≤ %.2f).\n', ...
        height(T), pval_t);
end

end % i_xct
