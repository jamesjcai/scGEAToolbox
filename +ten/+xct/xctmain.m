function [T] = xctmain(X_s, X_t, g, varargin)
%XCTMAIN  Lightweight spectral cell-cell interaction analysis.
%   Pure base-MATLAB implementation — no extra toolboxes required.
%   Uses Pearson co-expression as a GRN proxy instead of sc_pcnetpar,
%   making it fully self-contained.  For the full PCNet-backed version
%   see ten.sctenifoldxct.
%
%   T        = ten.xct.xctmain(X_s, X_t, g)
%   [T1, T2] = ten.xct.xctmain(X_s, X_t, g, 'twosided', true)
%   T        = ten.xct.xctmain(X_s, X_t, g, Name, Value)
%
%   Inputs:
%     X_s  - genes × source-cells raw count (or normalised) matrix
%     X_t  - genes × target-cells raw count (or normalised) matrix
%     g    - gene name list matching rows of X_s / X_t (string or cell)
%
%   Name-Value pairs:
%     'twosided'  - run both source→target AND target→source (default: true)
%     'n_dim'     - spectral embedding dimensions (default: 50)
%     'mu'        - cross-type correspondence weight (default: 0.9)
%     'corr_thr'  - |Pearson r| cut-off to keep GRN edges (default: 0.3)
%     'pval'      - p-value cut-off; 1.0 returns all pairs (default: 1.0)
%     'verbose'   - print progress (default: true)
%
%   Output:
%     T  - table with columns: ligand, receptor, dist, p_value
%          When twosided=true, T is a cell {T1, T2}.
%
%   Algorithm  (Path B — spectral approximation, adapted from ten.i_ma):
%     1. Log-normalise counts (library-size then log1p).
%     2. Build within-type co-expression W by thresholded Pearson correlation.
%     3. Build sparse L-R correspondence matrix W12 from built-in database.
%     4. Assemble block weight matrix W = [W_s, μW12; μW12', W_t] and its
%        graph Laplacian L = diag(|W|·1) − W  (identical to ten.i_ma).
%     5. Compute the n_dim smallest non-trivial eigenvectors of L via eigs.
%     6. Rank L-R candidate pairs by Euclidean distance in embedding space.
%     7. Nonparametric left-tail null test against sampled background pairs.
%
%   No dependencies beyond base MATLAB (R2019b+).
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742

% ── Parse inputs ─────────────────────────────────────────────────────────
p = inputParser;
addOptional(p, 'twosided', true,  @islogical);
addOptional(p, 'n_dim',    50,    @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',       0.9,   @(x) isnumeric(x) && x > 0);
addOptional(p, 'corr_thr', 0.3,   @(x) isnumeric(x) && x >= 0 && x <= 1);
addOptional(p, 'pval',     1.0,   @(x) isnumeric(x) && x >= 0);
addOptional(p, 'verbose',  true,  @islogical);
parse(p, varargin{:});

twosided = p.Results.twosided;
n_dim    = round(p.Results.n_dim);
mu_      = p.Results.mu;
corr_thr = p.Results.corr_thr;
pval_t   = p.Results.pval;
verbose  = p.Results.verbose;

X_s = double(full(X_s));
X_t = double(full(X_t));
g   = string(g(:));

ng = size(X_s, 1);
if size(X_t, 1) ~= ng || numel(g) ~= ng
    error('xctmain: X_s, X_t and g must all have the same number of rows.');
end

if verbose
    fprintf('[xctmain] Source: %d genes × %d cells\n', ng, size(X_s,2));
    fprintf('[xctmain] Target: %d genes × %d cells\n', ng, size(X_t,2));
end

% ── Load L-R database ─────────────────────────────────────────────────────
pw1    = fileparts(mfilename('fullpath'));
lr_mat = fullfile(pw1, '..', '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.mat');
lr_txt = fullfile(pw1, '..', '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.txt');

if exist(lr_mat, 'file')
    db     = load(lr_mat, 'ligand', 'receptor');
    lig_db = upper(string(db.ligand(:)));
    rec_db = upper(string(db.receptor(:)));
else
    Tlr    = readtable(lr_txt, 'FileType', 'text', 'Delimiter', '\t');
    lig_db = upper(string(Tlr{:,3}));
    rec_db = upper(string(Tlr{:,5}));
end
if verbose
    fprintf('[xctmain] L-R database: %d pairs.\n', numel(lig_db));
end

% ── Log-normalise (library-size then log1p) ───────────────────────────────
X_s = i_lognorm(X_s);
X_t = i_lognorm(X_t);

% ── Run alignment ─────────────────────────────────────────────────────────
T1 = i_align(X_s, X_t, g, lig_db, rec_db, n_dim, mu_, corr_thr, pval_t, verbose);

if twosided
    T2 = i_align(X_t, X_s, g, lig_db, rec_db, n_dim, mu_, corr_thr, pval_t, verbose);
    T  = {T1, T2};
else
    T = T1;
end

end % xctmain


%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────

function T = i_align(X_s, X_t, g, lig_db, rec_db, n_dim, mu_, corr_thr, pval_t, verbose)
%I_ALIGN  Core spectral alignment and null test for one direction.

ng   = size(X_s, 1);
g_up = upper(g);

% ── Step 1: within-type co-expression (base MATLAB Pearson, no toolbox) ──
if verbose, fprintf('[xctmain]   Computing co-expression matrices ...\n'); end
W11 = i_coexpr(X_s, corr_thr);   % sparse, ng × ng
W22 = i_coexpr(X_t, corr_thr);   % sparse, ng × ng

% ── Step 2: L-R correspondence matrix W12 (sparse, ng × ng) ───────────────
n_lr   = numel(lig_db);
li_buf = zeros(n_lr, 1);
ri_buf = zeros(n_lr, 1);
n_hit  = 0;
for k = 1:n_lr
    li = find(g_up == lig_db(k), 1);
    ri = find(g_up == rec_db(k), 1);
    if ~isempty(li) && ~isempty(ri)
        n_hit = n_hit + 1;
        li_buf(n_hit) = li;
        ri_buf(n_hit) = ri;
    end
end
li_idx = li_buf(1:n_hit);
ri_idx = ri_buf(1:n_hit);

if n_hit == 0
    warning('xctmain:noPairs', 'No L-R pairs found in gene list.');
    T = table(); return;
end
if verbose
    fprintf('[xctmain]   %d L-R pairs matched.\n', n_hit);
end
W12 = sparse(li_idx, ri_idx, ones(n_hit,1), ng, ng);

% ── Step 3: block weight matrix (follows ten.i_ma convention exactly) ─────
%   W = [W11+1,   mu*W12  ]
%       [mu*W12', W22+1   ]
%
%   Adding 1 to within-type matrices (as in i_ma.m) ensures the combined
%   graph is connected so the Laplacian has a unique trivial eigenvector.
W11_shift = W11 + speye(ng);   % diagonal +1 (sparse, avoids dense fill)
W22_shift = W22 + speye(ng);

w1_sum  = full(sum(abs(W11_shift(:))));
w2_sum  = full(sum(abs(W22_shift(:))));
w12_sum = full(sum(abs(W12(:))));
mu_scale = mu_ * (w1_sum + w2_sum) / (2 * max(w12_sum, eps));

W = [W11_shift,               mu_scale .* W12;
     mu_scale .* W12',        W22_shift      ];

% ── Step 4: sparse graph Laplacian ────────────────────────────────────────
d = full(sum(abs(W), 2));
L = spdiags(d, 0, 2*ng, 2*ng) - W;

% ── Step 5: spectral embedding ────────────────────────────────────────────
if verbose, fprintf('[xctmain]   Running eigs (n_dim=%d) ...\n', n_dim); end

n_ev = min(n_dim + 4, 2*ng - 2);
opts_eigs.isreal = true;
opts_eigs.issym  = true;
opts_eigs.tol    = 1e-6;
[V, D_ev] = eigs(L, n_ev, 'smallestreal', opts_eigs);
ev = real(diag(D_ev));
V  = real(V);
[ev, ord] = sort(ev, 'ascend');
V = V(:, ord);

% Discard trivial (near-zero) eigenvectors
keep = ev >= 1e-8;
V    = V(:, keep);

n_use = min(n_dim, size(V,2));
if n_use < n_dim && verbose
    warning('xctmain:dimReduced', ...
        'Only %d non-trivial eigenvectors found; n_dim reduced.', n_use);
end
V = V(:, 1:n_use);           % (2*ng) × n_use

P_s = V(1:ng,    :);         % source gene embeddings
P_t = V(ng+1:end,:);         % target gene embeddings

% ── Step 6: candidate L-R pair distances ─────────────────────────────────
cand_d = zeros(n_hit, 1);
for k = 1:n_hit
    d_k       = P_s(li_idx(k),:) - P_t(ri_idx(k),:);
    cand_d(k) = sqrt(d_k * d_k');
end

% ── Step 7: null distribution (sampled random non-L-R pairs) ─────────────
n_null  = max(50000, 100 * n_hit);
n_null  = min(n_null, ng^2 - n_hit);
rand_i  = randi(ng, n_null, 1);
rand_j  = randi(ng, n_null, 1);
diff_n  = P_s(rand_i,:) - P_t(rand_j,:);
null_d  = sqrt(sum(diff_n.^2, 2));

% Left-tail p-value: fraction of null ≤ candidate (small = strong signal)
p_vals = sum(null_d(:) <= cand_d(:)', 1)' ./ numel(null_d);

% ── Assemble output ───────────────────────────────────────────────────────
w12_vals = full(W12(sub2ind([ng,ng], li_idx, ri_idx)));
T = table(g(li_idx), g(ri_idx), cand_d, w12_vals, p_vals, ...
    'VariableNames', {'ligand','receptor','dist','correspondence','p_value'});
T = T(T.p_value <= pval_t, :);
T = sortrows(T, 'dist', 'ascend');

if verbose
    fprintf('[xctmain]   %d pairs returned (pval ≤ %.2f).\n', height(T), pval_t);
end

end % i_align


function W = i_coexpr(X, thr)
%I_COEXPR  Thresholded Pearson co-expression matrix (sparse).
%   Uses only base MATLAB — no Statistics Toolbox.
%   X : genes × cells.  Returns sparse ng × ng matrix.

% Log-scale already applied upstream; just z-score each gene across cells.
mu  = mean(X, 2);
sig = std(X, 0, 2);        % base MATLAB std (1/(n-1) normalisation)
sig(sig < eps) = 1;         % guard zero-variance genes
Xz  = (X - mu) ./ sig;     % genes × cells, each gene has unit std

nc = size(X, 2);
C  = (Xz * Xz') ./ (nc - 1);   % Pearson correlation matrix (ng × ng)
C  = max(min(C, 1), -1);        % numerical clamp

% Zero out weak and self correlations, then sparsify
C(abs(C) < thr) = 0;
C(1:size(C,1)+1:end) = 0;       % zero diagonal (self-correlation)
W = sparse(C);

end % i_coexpr


function X = i_lognorm(X)
%I_LOGNORM  Library-size normalisation followed by log1p.

col_sums = sum(X, 1);
col_sums(col_sums == 0) = 1;    % avoid /0 for empty cells
X = X ./ col_sums .* median(col_sums);
X = log1p(X);

end % i_lognorm
