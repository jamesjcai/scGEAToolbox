function [T, emb] = xctmain(X_s, X_t, g, varargin)
% XCTMAIN  Lightweight spectral cell-cell interaction analysis.
%   Builds within-type networks with PCNet (ten.i_xctgrn), configured to match
%   the published Python scTenifoldXct so the two are comparable.  For the
%   full-featured entry point, including the glyco channels, see
%   ten.sctenifoldxct.
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
%     'n_dim'     - spectral embedding dimensions (default: 3, as the reference)
%     'mu'        - cross-type correspondence weight (default: 1.0, as the reference)
%     'ncomp'     - PCNet principal components (default: 5, as in core.py)
%     'grn_q'     - GRN edge-filtering quantile; 0 disables (default: 0)
%     'pval'      - p-value cut-off; 1.0 returns all pairs (default: 1.0)
%     'verbose'   - print progress (default: true)
%
%   Output:
%     T  - table with columns: ligand, receptor, dist, p_value
%          When twosided=true, T is a cell {T1, T2}.
%
%   Algorithm  (Path B — spectral approximation, adapted from ten.i_ma):
%     1. Log-normalise counts (library-size then log1p).
%     2. Build within-type PCNet GRNs via ten.i_xctgrn.
%     3. Build sparse L-R correspondence matrix W12 from built-in database.
%     4. Assemble block weight matrix W = [W_s, μW12; μW12', W_t] and its
%        graph Laplacian L = diag(|W|·1) − W  (identical to ten.i_ma).
%     5. Compute the n_dim smallest non-trivial eigenvectors of L via eigs.
%     6. Rank L-R candidate pairs by Euclidean distance in embedding space.
%     7. Nonparametric left-tail null test against sampled background pairs.
%
%   Requires Statistics and Machine Learning Toolbox, via net.pcrnet.
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742

% ── Parse inputs ─────────────────────────────────────────────────────────
i_rejectcorrthr(varargin);

p = inputParser;
addOptional(p, 'twosided', true,  @islogical);
addOptional(p, 'n_dim',    3,     @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',       1.0,   @(x) isnumeric(x) && x > 0);
addOptional(p, 'ncomp',    5,     @(x) isnumeric(x) && x >= 2);
addOptional(p, 'grn_q',    0,     @(x) isnumeric(x) && x >= 0 && x < 1);
addOptional(p, 'pval',     1.0,   @(x) isnumeric(x) && x >= 0);
addOptional(p, 'verbose',  true,  @islogical);
addParameter(p, 'w12mode', "outer", @(x) ismember(string(x), ["outer","lr"]));
addParameter(p, 'alpha',   0.5,   @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'grnoffset', 1.0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'useparallel', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});

twosided = p.Results.twosided;
n_dim    = round(p.Results.n_dim);
mu_      = p.Results.mu;
ncomp    = round(p.Results.ncomp);
grn_q    = p.Results.grn_q;
pval_t   = p.Results.pval;
verbose  = p.Results.verbose;

cfg = struct('n_dim', n_dim, 'mu', mu_, 'ncomp', ncomp, 'grn_q', grn_q, ...
    'pval', pval_t, 'verbose', verbose, ...
    'w12mode', string(p.Results.w12mode), 'alpha', p.Results.alpha, ...
    'grnoffset', p.Results.grnoffset, ...
    'useparallel', logical(p.Results.useparallel));

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
[T1, E1] = i_align(X_s, X_t, g, lig_db, rec_db, cfg);

if twosided
    [T2, E2] = i_align(X_t, X_s, g, lig_db, rec_db, cfg);
    T  = {T1, T2};
    emb = {E1, E2};
else
    T = T1;
    emb = E1;
end

end % xctmain


%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────

function [T, emb] = i_align(X_s, X_t, g, lig_db, rec_db, cfg)
% I_ALIGN  Core spectral alignment and null test for one direction.

n_dim   = cfg.n_dim;
mu_     = cfg.mu;
ncomp   = cfg.ncomp;
grn_q   = cfg.grn_q;
pval_t  = cfg.pval;
verbose = cfg.verbose;
w12mode = cfg.w12mode;
alpha_  = cfg.alpha;
grnoffset = cfg.grnoffset;
useparallel = cfg.useparallel;

ng   = size(X_s, 1);
g_up = upper(g);

% ── Step 1: within-type PCNet GRNs ───────────────────────────────────────
if verbose, fprintf('[xctmain]   Building PCNet GRNs ...\n'); end
W11 = ten.i_xctgrn(X_s, ncomp, grn_q, verbose, useparallel);   % sparse, ng × ng
W22 = ten.i_xctgrn(X_t, ncomp, grn_q, verbose, useparallel);   % sparse, ng × ng

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

emb = struct('P_s', [], 'P_t', [], 'li_idx', li_idx, 'ri_idx', ri_idx);
if n_hit == 0
    warning('xctmain:noPairs', 'No L-R pairs found in gene list.');
    T = table(); return;
end
if verbose
    fprintf('[xctmain]   %d L-R pairs matched.\n', n_hit);
end
% Correspondence block. "outer" reproduces core.py:_build_w with its default
% query_DB=None, in which the L-R database does not enter the alignment at all
% and is applied only when candidates are ranked below.
W12 = ten.i_xctw12(X_s, X_t, ng, mode=w12mode, alpha=alpha_, ...
    li_idx=li_idx, ri_idx=ri_idx, lognorm=false);

% ── Step 3: block weight matrix ───────────────────────────────────────────
%   W = [W11+offset,  mu*W12    ]
%       [mu*W12',     W22+offset]
%
%   The offset is added to EVERY element, as core.py:415 and ten.i_ma do. This
%   file previously wrote W11 + speye(ng), a diagonal-only shift, which is a
%   no-op for the Laplacian and left the graph disconnected. See TEN.I_XCTBLOCK.
W = ten.i_xctblock(W11, W22, W12, mu_, grnoffset);

% ── Steps 4-5: graph Laplacian and spectral embedding ─────────────────────
if verbose, fprintf('[xctmain]   Running eigs (n_dim=%d) ...\n', n_dim); end

[P_s, P_t] = ten.i_laplacianembed(W, ng, n_dim, verbose);

% Returned so ten.xct.xctmain2 can build its chi-square statistic over all gene
% pairs, which the candidate-only output table cannot support.
emb.P_s = P_s;
emb.P_t = P_t;

% ── Step 6: candidate L-R pair distances ─────────────────────────────────
cand_d = zeros(n_hit, 1);
for k = 1:n_hit
    d_k       = P_s(li_idx(k),:) - P_t(ri_idx(k),:);
    cand_d(k) = sqrt(d_k * d_k');
end

% ── Step 7: left-tail null test ──────────────────────────────────────────
% Exhaustive over non-candidate pairs, as stat.py::null_test. Small p-value
% means the candidate distance is unusually small. See TEN.I_XCTNULL for why
% this replaced a sampled null.
p_vals = ten.i_xctnull(P_s, P_t, li_idx, ri_idx, W12);

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


function i_rejectcorrthr(args)
% I_REJECTCORRTHR  Explain the removal of the Pearson GRN option.
%   The 'corr_thr' argument selected a cut-off for a Pearson co-expression
%   proxy that no longer exists, so silently ignoring it would leave callers
%   believing they had configured a network they had not.

isname = cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, 'corr_thr'), args);
if any(isname)
    error("TEN:XCT:XCTMAIN:CorrThrRemoved", ...
        "'corr_thr' is no longer supported. Within-type networks are now " + ...
        "built with PCNet rather than thresholded Pearson co-expression, " + ...
        "to match the published Python implementation. Use 'ncomp' to set " + ...
        "the number of principal components, or 'grn_q' to filter edges.");
end

end % i_rejectcorrthr


function X = i_lognorm(X)
% I_LOGNORM  Library-size normalisation followed by log1p.

col_sums = sum(X, 1);
col_sums(col_sums == 0) = 1;    % avoid /0 for empty cells
X = X ./ col_sums .* median(col_sums);
X = log1p(X);

end % i_lognorm
