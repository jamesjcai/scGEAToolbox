function [T] = sctenifoldxct2(sce1, sce2, celltype1, celltype2, twosided, varargin)
%SCTENIFOLDXCT2  Differential cell-cell interaction across two samples.
%   Native MATLAB implementation (Path B, PCNet GRNs).  No Python required.
%
%   T = ten.sctenifoldxct2(SCE1, SCE2, CT1, CT2)
%   T = ten.sctenifoldxct2(SCE1, SCE2, CT1, CT2, TWOSIDED)
%   T = ten.sctenifoldxct2(SCE1, SCE2, CT1, CT2, TWOSIDED, Name, Value)
%
%   Identifies L-R pairs whose interaction strength differs significantly
%   between two conditions (SCE1 vs SCE2) for the same pair of cell types.
%   For each direction, spectral manifold alignment is run independently on
%   each sample; then fold-change of embedding distances is tested against a
%   null distribution of random (non-L-R) gene pair FCs.
%
%   Inputs:
%     SCE1, SCE2  - SingleCellExperiment objects for the two conditions
%     CT1         - source cell type label (string)
%     CT2         - target cell type label (string)
%     TWOSIDED    - logical, run CT1→CT2 AND CT2→CT1 (default: true)
%
%   Name-Value pairs:
%     'n_dim'   - spectral embedding dimensions (default: 50)
%     'mu'      - cross-type correspondence weight (default: 0.9)
%     'verbose' - print progress messages (default: true)
%
%   Output:
%     T  - table: ligand, receptor, dist_s1, dist_s2, FC, p_value
%          FC < 1 : interaction stronger in SCE1 (gained in condition 1)
%          FC > 1 : interaction stronger in SCE2 (gained in condition 2)
%          p_value: two-tailed empirical significance of |log2(FC)|
%          When twosided=true, T is a cell array {T_ct1_ct2, T_ct2_ct1}.
%
%   This is a MATLAB replacement for +run/py_scTenifoldXct2.
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742

if nargin < 5, twosided = true; end

p = inputParser;
addOptional(p, 'n_dim',   50,   @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',      0.9,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'verbose', true, @islogical);
parse(p, varargin{:});

n_dim   = round(p.Results.n_dim);
mu_     = p.Results.mu;
verbose = p.Results.verbose;

% ── 1. Prepare expression matrices ───────────────────────────────────────
[X_s1, X_t1, X_s2, X_t2, g] = i_prepare(sce1, sce2, celltype1, celltype2, verbose);

% ── 2. Load L-R database ─────────────────────────────────────────────────
[lig_db, rec_db] = i_loadlr();
if verbose
    fprintf('[sctenifoldxct2] L-R database: %d pairs loaded.\n', numel(lig_db));
end

% ── 3. Match L-R pairs to gene list ──────────────────────────────────────
g_up = upper(string(g(:)));
[li_idx, ri_idx] = i_matchlr(g_up, lig_db, rec_db);
n_valid = numel(li_idx);
if n_valid == 0
    warning('sctenifoldxct2:noPairs', 'No L-R pairs found in gene list.');
    T = table(); return;
end
if verbose
    fprintf('[sctenifoldxct2] %d L-R pairs matched in data.\n', n_valid);
end

% ── 4. Build GRNs for both samples ───────────────────────────────────────
if verbose, fprintf('[sctenifoldxct2] Building GRNs for sample 1...\n'); end
A_s1 = i_buildgrn(X_s1);
A_t1 = i_buildgrn(X_t1);
if verbose, fprintf('[sctenifoldxct2] Building GRNs for sample 2...\n'); end
A_s2 = i_buildgrn(X_s2);
A_t2 = i_buildgrn(X_t2);

% ── 5. Run differential alignment ────────────────────────────────────────
if verbose
    fprintf('[sctenifoldxct2] Differential alignment %s → %s ...\n', ...
        celltype1, celltype2);
end
T1 = i_xct2dir(g, A_s1, A_t1, A_s2, A_t2, li_idx, ri_idx, n_dim, mu_, verbose);

if twosided
    if verbose
        fprintf('[sctenifoldxct2] Differential alignment %s → %s ...\n', ...
            celltype2, celltype1);
    end
    T2 = i_xct2dir(g, A_t1, A_s1, A_t2, A_s2, li_idx, ri_idx, n_dim, mu_, verbose);
    T = {T1, T2};
else
    T = T1;
end

end % main function


%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────

function [X_s1, X_t1, X_s2, X_t2, g] = i_prepare(sce1, sce2, ct1, ct2, verbose)
    function [Xs, Xt] = from_sce(sce, c1, c2)
        X = sce.X;
        if issparse(X), X = full(X); end
        X  = single(X);
        Xs = X(:, sce.c_cell_type_tx == c1);
        Xt = X(:, sce.c_cell_type_tx == c2);
    end
    [X_s1, X_t1] = from_sce(sce1, ct1, ct2);
    [X_s2, X_t2] = from_sce(sce2, ct1, ct2);
    g = sce1.g;
    if verbose
        fprintf('[sctenifoldxct2] S1  %s: %d cells  |  %s: %d cells\n', ...
            ct1, size(X_s1,2), ct2, size(X_t1,2));
        fprintf('[sctenifoldxct2] S2  %s: %d cells  |  %s: %d cells\n', ...
            ct1, size(X_s2,2), ct2, size(X_t2,2));
    end
end

function A = i_buildgrn(X)
    A = sc_pcnetpar(X);
    A = A ./ max(abs(A(:)));
    A = ten.e_filtadjc(A, 0.75, false);
end

function [lig_db, rec_db] = i_loadlr()
    pw1    = fileparts(mfilename('fullpath'));
    lr_mat = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.mat');
    lr_txt = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.txt');
    if exist(lr_mat, 'file')
        db     = load(lr_mat, 'ligand', 'receptor');
        lig_db = upper(string(db.ligand(:)));
        rec_db = upper(string(db.receptor(:)));
    else
        T_lr   = readtable(lr_txt, 'FileType', 'text', 'Delimiter', '\t');
        lig_db = upper(string(T_lr{:,3}));
        rec_db = upper(string(T_lr{:,5}));
    end
end

function [li_idx, ri_idx] = i_matchlr(g_up, lig_db, rec_db)
    n_lr  = numel(lig_db);
    li    = zeros(n_lr, 1);
    ri    = zeros(n_lr, 1);
    n     = 0;
    for k = 1:n_lr
        l = find(g_up == lig_db(k), 1);
        r = find(g_up == rec_db(k), 1);
        if ~isempty(l) && ~isempty(r)
            n = n + 1;
            li(n) = l;
            ri(n) = r;
        end
    end
    li_idx = li(1:n);
    ri_idx = ri(1:n);
end

function Td = i_xct2dir(g, A_s1, A_t1, A_s2, A_t2, li_idx, ri_idx, n_dim, mu_, verbose)
%I_XCT2DIR  Differential alignment for one direction (source → target).

    ng      = size(A_s1, 1);
    n_valid = numel(li_idx);

    % Spectral embeddings for each sample
    [P_s1, P_t1] = i_embed(A_s1, A_t1, li_idx, ri_idx, ng, n_dim, mu_);
    [P_s2, P_t2] = i_embed(A_s2, A_t2, li_idx, ri_idx, ng, n_dim, mu_);

    % L-R embedding distances in each sample
    d1 = zeros(n_valid, 1, 'single');
    d2 = zeros(n_valid, 1, 'single');
    for k = 1:n_valid
        dv1   = P_s1(li_idx(k),:) - P_t1(ri_idx(k),:);
        dv2   = P_s2(li_idx(k),:) - P_t2(ri_idx(k),:);
        d1(k) = sqrt(sum(dv1.^2));
        d2(k) = sqrt(sum(dv2.^2));
    end

    % Null distribution: |log2(FC)| for random non-L-R gene pairs
    n_null    = max(50000, 100 * n_valid);
    rng_state = rng;
    rand_i    = randi(ng, n_null, 1);
    rand_j    = randi(ng, n_null, 1);
    rng(rng_state);
    null_d1   = sqrt(sum((single(P_s1(rand_i,:)) - single(P_t1(rand_j,:))).^2, 2));
    null_d2   = sqrt(sum((single(P_s2(rand_i,:)) - single(P_t2(rand_j,:))).^2, 2));
    null_lfc  = log2(null_d1 ./ max(null_d2, 1e-8));

    % Two-tailed differential p-value
    FC     = d1 ./ max(d2, 1e-8);
    lfc    = log2(FC);
    p_vals = sum(abs(null_lfc) >= abs(lfc)', 1)' ./ n_null;

    % Output table
    lig_out = g(li_idx);
    rec_out = g(ri_idx);
    Td = table(lig_out, rec_out, double(d1), double(d2), double(FC), p_vals, ...
        'VariableNames', {'ligand','receptor','dist_s1','dist_s2','FC','p_value'});
    Td = sortrows(Td, 'p_value', 'ascend');

    if verbose
        fprintf('[sctenifoldxct2]   %d L-R pairs; %d with p < 0.05.\n', ...
            height(Td), sum(Td.p_value < 0.05));
    end
end

function [P_s, P_t] = i_embed(A_s, A_t, li_idx, ri_idx, ng, n_dim, mu_)
%I_EMBED  Spectral manifold alignment (identical kernel to ten.sctenifoldxct).

    n_valid = numel(li_idx);
    W11     = sparse(0.5*(A_s + A_s'));
    W22     = sparse(0.5*(A_t + A_t'));
    W12     = sparse(li_idx, ri_idx, ones(n_valid,1), ng, ng);

    w11_s = sum(abs(W11(:)));
    w22_s = sum(abs(W22(:)));
    w12_s = sum(abs(W12(:)));
    if w12_s == 0
        mu_scale = 0;
    else
        mu_scale = mu_ * (w11_s + w22_s) / (2 * w12_s);
    end

    W = [W11,               mu_scale .* W12;
         mu_scale .* W12',  W22            ];
    d = full(sum(abs(W), 2));
    L = spdiags(d, 0, 2*ng, 2*ng) - W;

    n_ev = min(n_dim + 4, 2*ng - 2);
    opts.isreal = true;
    opts.issym  = true;
    opts.tol    = 1e-6;
    [V, D_ev] = eigs(L, n_ev, 'smallestreal', opts);
    ev = real(diag(D_ev));
    V  = real(V);
    [ev, ord] = sort(ev, 'ascend');
    V  = V(:, ord);
    keep = ev >= 1e-8;
    V    = V(:, keep);
    n_use = min(n_dim, size(V, 2));
    V     = V(:, 1:n_use);

    P_s = V(1:ng,    :);
    P_t = V(ng+1:end,:);
end
