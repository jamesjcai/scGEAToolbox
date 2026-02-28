function [T] = xctmain_nn(X_s, X_t, g, varargin)
%XCTMAIN_NN  Neural-network manifold alignment for cell-cell interaction.
%   Path A: faithful MATLAB translation of the published scTenifoldXct
%   training loop.  Requires Deep Learning Toolbox (R2022b+).
%
%   For a toolbox-free alternative see ten.xct.xctmain (Path B, spectral).
%
%   T        = ten.xct.xctmain_nn(X_s, X_t, g)
%   [T1, T2] = ten.xct.xctmain_nn(X_s, X_t, g, 'twosided', true)
%   T        = ten.xct.xctmain_nn(X_s, X_t, g, Name, Value)
%
%   Inputs:
%     X_s  - genes × source-cells expression matrix
%     X_t  - genes × target-cells expression matrix
%     g    - gene name list matching rows of X_s / X_t (string or cell)
%
%   Name-Value pairs:
%     'twosided'  - run both directions (default: true)
%     'n_dim'     - embedding dimension (default: 50)
%     'mu'        - cross-type correspondence weight (default: 0.9)
%     'corr_thr'  - |Pearson r| cut-off for co-expression edges (default: 0.3)
%     'n_steps'   - Adam training iterations (default: 1000)
%     'lr'        - Adam learning rate (default: 0.01)
%     'pval'      - p-value cut-off; 1.0 returns all pairs (default: 1.0)
%     'verbose'   - print progress (default: true)
%
%   Output:
%     T  - table: ligand, receptor, dist, p_value
%          When twosided=true, T is a cell {T1, T2}.
%
%   Algorithm (Path A — neural network, faithful to Python scTenifoldXct):
%
%     Two 3-layer feedforward networks f_s and f_t (one per cell type) map
%     each gene's expression profile across cells to an n_dim-dimensional
%     embedding.  Architecture for each network:
%
%       input (n_cells)
%         → FC(4·n_h) + sigmoid          n_h = floor(sqrt(n_cells · n_dim))
%         → FC(n_h)   + sigmoid
%         → FC(n_dim)
%
%     Training minimises the graph-Laplacian loss:
%
%       L = trace(P' · L_W · P) / 3000
%
%     where P = U·V'  (Stiefel retraction via economy SVD of stacked outputs)
%     and L_W is the unnormalised Laplacian of the block weight matrix
%
%       W = [ GRN_source,   μ·W_LR  ]
%           [ μ·W_LR',   GRN_target ]
%
%     Ligand-receptor pairs in W_LR pull the corresponding gene embeddings
%     together; the GRN structure keeps biologically related genes adjacent.
%     Gradients flow through dlsvd; both networks are updated with Adam.
%
%   This is mathematically equivalent to the Python implementation.
%   The GRN proxy here is thresholded Pearson co-expression instead of
%   sc_pcnetpar; for the PCNet-backed version see ten.sctenifoldxct.
%
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742

% ── Toolbox guard ─────────────────────────────────────────────────────────
if ~(exist('dlarray', 'builtin') || exist('dlarray', 'file'))
    error(['xctmain_nn requires the MATLAB Deep Learning Toolbox (R2022b+). ' ...
           'Use ten.xct.xctmain for a toolbox-free spectral alternative.']);
end

% ── Parse inputs ─────────────────────────────────────────────────────────
p = inputParser;
addOptional(p, 'twosided', true,  @islogical);
addOptional(p, 'n_dim',    50,    @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',       0.9,   @(x) isnumeric(x) && x > 0);
addOptional(p, 'corr_thr', 0.3,   @(x) isnumeric(x) && x >= 0 && x <= 1);
addOptional(p, 'n_steps',  1000,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'lr',       0.01,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'pval',     1.0,   @(x) isnumeric(x) && x >= 0);
addOptional(p, 'verbose',  true,  @islogical);
parse(p, varargin{:});

twosided = p.Results.twosided;
n_dim    = round(p.Results.n_dim);
mu_      = p.Results.mu;
corr_thr = p.Results.corr_thr;
n_steps  = round(p.Results.n_steps);
lr       = p.Results.lr;
pval_t   = p.Results.pval;
verbose  = p.Results.verbose;

X_s = double(full(X_s));
X_t = double(full(X_t));
g   = string(g(:));

ng = size(X_s, 1);
if size(X_t, 1) ~= ng || numel(g) ~= ng
    error('xctmain_nn: X_s, X_t and g must all have the same number of rows.');
end
if ng > 3000
    warning('xctmain_nn:largeGeneSet', ...
        ['%d genes → L matrix will be %d × %d (%.0f MB single). ' ...
         'Consider subsetting to top HVGs to reduce memory.'], ...
        ng, 2*ng, 2*ng, (2*ng)^2*4/1e6);
end

if verbose
    fprintf('[xctmain_nn] Source: %d genes × %d cells\n', ng, size(X_s,2));
    fprintf('[xctmain_nn] Target: %d genes × %d cells\n', ng, size(X_t,2));
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
    fprintf('[xctmain_nn] L-R database: %d pairs.\n', numel(lig_db));
end

% ── Log-normalise (library-size then log1p) ───────────────────────────────
X_s = i_lognorm(X_s);
X_t = i_lognorm(X_t);

% ── Run ───────────────────────────────────────────────────────────────────
T1 = i_align_nn(X_s, X_t, g, lig_db, rec_db, ...
                n_dim, mu_, corr_thr, n_steps, lr, pval_t, verbose);

if twosided
    T2 = i_align_nn(X_t, X_s, g, lig_db, rec_db, ...
                    n_dim, mu_, corr_thr, n_steps, lr, pval_t, verbose);
    T  = {T1, T2};
else
    T = T1;
end

end % xctmain_nn


%% ── CORE ALIGNMENT ───────────────────────────────────────────────────────

function T = i_align_nn(X_s, X_t, g, lig_db, rec_db, ...
                        n_dim, mu_, corr_thr, n_steps, lr, pval_t, verbose)
%I_ALIGN_NN  Neural-network manifold alignment for one direction.

ng        = size(X_s, 1);
n_cells_s = size(X_s, 2);
n_cells_t = size(X_t, 2);
g_up      = upper(g);

% ── 1. Pearson co-expression (within-type GRN proxy) ─────────────────────
if verbose, fprintf('[xctmain_nn]   Building co-expression matrices ...\n'); end
W11 = i_coexpr(X_s, corr_thr);   % sparse ng × ng
W22 = i_coexpr(X_t, corr_thr);   % sparse ng × ng

% ── 2. L-R correspondence matrix (sparse ng × ng) ────────────────────────
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
    warning('xctmain_nn:noPairs', 'No L-R pairs found. Returning empty table.');
    T = table(); return;
end
if verbose, fprintf('[xctmain_nn]   %d L-R pairs matched.\n', n_hit); end

W12 = sparse(li_idx, ri_idx, ones(n_hit,1), ng, ng);

% ── 3. Block weight matrix W and graph Laplacian L ───────────────────────
%   Follows Python scTenifoldXct exactly — no +1 shift (raw GRN + L-R links).
w1_sum  = sum(abs(W11(:)));
w2_sum  = sum(abs(W22(:)));
w12_sum = sum(abs(W12(:)));
mu_scale = mu_ * (w1_sum + w2_sum) / (2 * max(w12_sum, eps));

W = [W11,               mu_scale .* W12;
     mu_scale .* W12',  W22            ];    % (2*ng) × (2*ng) sparse

% Unnormalised graph Laplacian: L = diag(|W|·1) − W
d = full(sum(abs(W), 2));
L_dense = single(full(spdiags(d, 0, 2*ng, 2*ng) - W));  % dense single

% ── 4. Initialise networks ────────────────────────────────────────────────
%   Architecture: D_in → FC(4·n_h,σ) → FC(n_h,σ) → FC(n_dim)
%   D_in = n_cells of that cell type; n_h = floor(sqrt(n_cells · n_dim))
rng(0);  % reproducibility
ps = i_init_params(n_cells_s, n_dim);   % source network
pt = i_init_params(n_cells_t, n_dim);   % target network

% Convert data to single-precision dlarray (unlabelled = custom training)
Xs_dl = dlarray(single(X_s));   % ng × n_cells_s
Xt_dl = dlarray(single(X_t));   % ng × n_cells_t
L_dl  = dlarray(L_dense);       % (2*ng) × (2*ng)  — constant w.r.t. params

% ── 5. Training loop (Adam, Stiefel retraction via SVD each step) ─────────
%   Mirrors Python ManifoldAlignmentNet.train() exactly:
%     outputs = cat([net_s(X_s); net_t(X_t)])
%     P = U·V'  from SVD(outputs)           ← Stiefel retraction
%     loss = trace(P'·L·P) / 3000
%     backprop + Adam step
if verbose
    fprintf('[xctmain_nn]   Training (%d steps, lr=%.4f) ...\n', n_steps, lr);
end

[avg_s,   avgSq_s]   = deal([]);   % Adam trailing-average state (auto-init)
[avg_t,   avgSq_t]   = deal([]);

for step = 1:n_steps
    [loss_val, grads] = dlfeval(@i_grad_fn, ps, pt, Xs_dl, Xt_dl, L_dl);

    grads_s = grads{1};
    grads_t = grads{2};

    [ps, avg_s, avgSq_s] = adamupdate(ps, grads_s, avg_s, avgSq_s, step, lr);
    [pt, avg_t, avgSq_t] = adamupdate(pt, grads_t, avg_t, avgSq_t, step, lr);

    if verbose && (step == 1 || mod(step, 100) == 0)
        fprintf('[xctmain_nn]     step %4d / %d   loss = %.6f\n', ...
            step, n_steps, extractdata(loss_val));
    end
end

% ── 6. Extract final embeddings ───────────────────────────────────────────
out_s   = extractdata(i_net_fwd(Xs_dl, ps));   % ng × n_dim
out_t   = extractdata(i_net_fwd(Xt_dl, pt));   % ng × n_dim
outputs = [out_s; out_t];                       % (2*ng) × n_dim

% Stiefel projection (same retraction as during training)
[U, ~, V] = svd(outputs, 'econ');
P   = U * V';      % (2*ng) × n_dim
P_s = P(1:ng, :);
P_t = P(ng+1:end, :);

% ── 7. Candidate L-R distances ────────────────────────────────────────────
cand_d = zeros(n_hit, 1);
for k = 1:n_hit
    d_k      = P_s(li_idx(k),:) - P_t(ri_idx(k),:);
    cand_d(k) = sqrt(d_k * d_k');
end

% ── 8. Null distribution (sampled random non-L-R pairs) ──────────────────
n_null = max(50000, 100 * n_hit);
n_null = min(n_null, ng^2 - n_hit);
rand_i = randi(ng, n_null, 1);
rand_j = randi(ng, n_null, 1);
diff_n = P_s(rand_i,:) - P_t(rand_j,:);
null_d = sqrt(sum(diff_n.^2, 2));

% Left-tail p-value: fraction of null ≤ candidate (small = strong signal)
p_vals = sum(null_d(:) <= cand_d(:)', 1)' ./ numel(null_d);

% ── 9. Assemble output ────────────────────────────────────────────────────
w12_vals = full(W12(sub2ind([ng,ng], li_idx, ri_idx)));
T = table(g(li_idx), g(ri_idx), cand_d, w12_vals, p_vals, ...
    'VariableNames', {'ligand','receptor','dist','correspondence','p_value'});
T = T(T.p_value <= pval_t, :);
T = sortrows(T, 'dist', 'ascend');

if verbose
    fprintf('[xctmain_nn]   %d pairs returned (pval ≤ %.2f).\n', ...
        height(T), pval_t);
end

end % i_align_nn


%% ── NEURAL NETWORK PRIMITIVES ────────────────────────────────────────────

function [loss, grads] = i_grad_fn(ps, pt, Xs_dl, Xt_dl, L_dl)
%I_GRAD_FN  Loss and gradients — called inside dlfeval.
%
%   loss = trace(P' · L · P) / 3000
%   where P = U·V'  (SVD Stiefel retraction of stacked net outputs).
%   Gradients flow back through dlsvd to both network parameter structs.

out_s   = i_net_fwd(Xs_dl, ps);   % ng_s × n_dim  (dlarray)
out_t   = i_net_fwd(Xt_dl, pt);   % ng_t × n_dim  (dlarray)
outputs = [out_s; out_t];          % (ng_s+ng_t) × n_dim

% Stiefel retraction: P = U·V'  (economy SVD)
% dlsvd supports automatic differentiation (R2022b+)
[U, ~, V] = dlsvd(outputs, 'econ');
P = U * V';                        % (ng_s+ng_t) × n_dim  (dlarray)

% Laplacian loss  trace(P'·L·P)/3000  =  sum(P .* (L·P)) / 3000
% L_dl is a constant dlarray — gradients accumulate only through P
loss = sum(P .* (L_dl * P), 'all') / 3000;

% Gradients w.r.t. both parameter structs simultaneously
grads = dlgradient(loss, {ps, pt});

end % i_grad_fn


function out = i_net_fwd(X, p)
%I_NET_FWD  Forward pass through one 3-layer network.
%
%   X   : ng × n_cells  (each row = one gene's expression profile)
%   p   : struct with fields W1,b1,W2,b2,W3,b3  (dlarray)
%   out : ng × n_dim    (each row = one gene's embedding)
%
%   Layer sizes (set at init): n_cells → 4·n_h → n_h → n_dim
%   Activations: sigmoid on hidden layers, linear output.

h1  = sigmoid(X  * p.W1 + p.b1);   % ng × H1
h2  = sigmoid(h1 * p.W2 + p.b2);   % ng × H2
out =          h2 * p.W3 + p.b3;   % ng × n_dim

end % i_net_fwd


function p = i_init_params(n_cells, n_dim)
%I_INIT_PARAMS  Xavier-initialised parameter struct for one network.
%   n_h = floor(sqrt(n_cells * n_dim))  (geometric mean, as in Python)
%   H1  = 4 * n_h,  H2 = n_h

n_h  = max(floor(sqrt(n_cells * n_dim)), 1);
H1   = 4 * n_h;
H2   = n_h;

% Xavier (Glorot) uniform initialisation  ±sqrt(6/(fan_in+fan_out))
p.W1 = dlarray(i_xavier(n_cells, H1,    'single'));
p.b1 = dlarray(zeros(1, H1,             'single'));
p.W2 = dlarray(i_xavier(H1,    H2,      'single'));
p.b2 = dlarray(zeros(1, H2,             'single'));
p.W3 = dlarray(i_xavier(H2,    n_dim,   'single'));
p.b3 = dlarray(zeros(1, n_dim,          'single'));

end % i_init_params


function W = i_xavier(fan_in, fan_out, prec)
%I_XAVIER  Xavier uniform weight matrix.
lim = sqrt(6 / (fan_in + fan_out));
W   = (rand(fan_in, fan_out, prec) * 2 - 1) .* lim;
end % i_xavier


%% ── SHARED HELPERS (duplicated from xctmain for self-containment) ────────

function W = i_coexpr(X, thr)
%I_COEXPR  Thresholded Pearson co-expression (sparse, base MATLAB only).
mu  = mean(X, 2);
sig = std(X, 0, 2);
sig(sig < eps) = 1;
Xz  = (X - mu) ./ sig;
C   = (Xz * Xz') ./ max(size(X,2) - 1, 1);
C   = max(min(C, 1), -1);
C(abs(C) < thr) = 0;
C(1:size(C,1)+1:end) = 0;   % zero diagonal
W   = sparse(C);
end % i_coexpr


function X = i_lognorm(X)
%I_LOGNORM  Library-size normalisation + log1p.
cs = sum(X, 1);
cs(cs == 0) = 1;
X  = log1p(X ./ cs .* median(cs));
end % i_lognorm
