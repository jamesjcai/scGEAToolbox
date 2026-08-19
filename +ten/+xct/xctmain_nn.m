function [T, emb] = xctmain_nn(X_s, X_t, g, varargin)
% XCTMAIN_NN  Neural-network manifold alignment for cell-cell interaction.
%   Path A: MATLAB translation of the published scTenifoldXct training loop
%   (see Provenance below).  Requires Deep Learning Toolbox (R2022b+).
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
%     'n_dim'     - embedding dimension (default: 3, as the reference)
%     'mu'        - cross-type correspondence weight (default: 1.0, as the reference)
%     'ncomp'     - PCNet principal components (default: 5, as in core.py)
%     'grn_q'     - GRN edge-filtering quantile; 0 disables (default: 0)
%     'n_steps'   - Adam training iterations (default: 1000)
%     'lr'        - Adam learning rate (default: 0.01)
%     'pval'      - p-value cut-off; 1.0 returns all pairs (default: 1.0)
%     'verbose'   - print progress (default: true)
%
%   Output:
%     T  - table: ligand, receptor, dist, p_value
%          When twosided=true, T is a cell {T1, T2}.
%
%   Algorithm (Path A — neural network, following Python scTenifoldXct):
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
%     where P is the Stiefel retraction of the stacked network outputs and
%     L_W is the unnormalised Laplacian of the block weight matrix
%
%       W = [ GRN_source,   μ·W_LR  ]
%           [ μ·W_LR',   GRN_target ]
%
%     Ligand-receptor pairs in W_LR pull the corresponding gene embeddings
%     together; the GRN structure keeps biologically related genes adjacent.
%     Both networks are updated with Adam.  Each step retracts the stacked
%     network outputs onto the Stiefel manifold with an exact economy SVD, so
%     P'*P = I holds throughout training, not only at the end.  ten.i_nnembed
%     owns the training loop and reports the orthogonality error of what it
%     returns.
%
%   Provenance: a reimplementation of the published Python scTenifoldXct,
%   written against that package rather than ported from vendored source
%   (external/py_scTenifoldXct holds only driver scripts that import it).
%   Checked line-by-line against nn.py, stiefel.py and core.py of
%   cailab-tamu/scTenifoldXct.
%
%   Verified identical: the 3-layer sigmoid/sigmoid/linear architecture; the
%   layer widths, since scipy gmean([n_cells, n_dim]) is sqrt(n_cells*n_dim)
%   and their a=4 is the 4*n_h here; torch.nn.Linear's default initialisation;
%   Adam at lr 0.01 for 1000 steps with matching betas and epsilon; the loss
%   trace(P'*L*P)/3000, including that unexplained /3000; the Laplacian degree
%   convention; the Riemannian gradient and its backpropagation through the
%   retraction; and the returned embedding, which the reference takes from
%   inside the final iteration rather than recomputing afterwards.
%
%   One deviation remains, and it is unavoidable: the RNG streams differ, so
%   the two cannot produce bitwise-identical weights from the same seed. Every
%   comparison is therefore distributional, read against the reference's own
%   seed-to-seed variation as the ceiling.
%
%   Note that MATLAB does not apply core.py's +1 GRN offset. That is a
%   caller-level choice, not part of this training loop, but it means W here is
%   signed where the reference's is not; see the DEGREE CONVENTION section of
%   TEN.I_NNEMBED for what that implies.
%
%   The GRN is PCNet via ten.i_xctgrn, configured to match core.py's
%   make_pcNet(nComp=5) with no edge filtering. Note that ten.sctenifoldxct
%   does NOT currently use those settings - it builds its networks with
%   ncomp=3 and a 0.75 quantile filter - so results from the two entry points
%   are not directly comparable.
%
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742

% ── Toolbox guard ─────────────────────────────────────────────────────────
if ~(exist('dlarray', 'builtin') || exist('dlarray', 'file'))
    error(['xctmain_nn requires the MATLAB Deep Learning Toolbox (R2022b+). ' ...
           'Use ten.xct.xctmain for a toolbox-free spectral alternative.']);
end

% ── Parse inputs ─────────────────────────────────────────────────────────
i_rejectcorrthr(varargin);

p = inputParser;
addParameter(p, 'twosided', true,  @islogical);
addParameter(p, 'n_dim',    3,     @(x) isnumeric(x) && x > 0);
addParameter(p, 'mu',       1.0,   @(x) isnumeric(x) && x > 0);
addParameter(p, 'ncomp',    5,     @(x) isnumeric(x) && x >= 2);
addParameter(p, 'grn_q',    0,     @(x) isnumeric(x) && x >= 0 && x < 1);
addParameter(p, 'n_steps',  1000,  @(x) isnumeric(x) && x > 0);
addParameter(p, 'lr',       0.01,  @(x) isnumeric(x) && x > 0);
addParameter(p, 'pval',     1.0,   @(x) isnumeric(x) && x >= 0);
addParameter(p, 'verbose',  true,  @islogical);
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
n_steps  = round(p.Results.n_steps);
lr       = p.Results.lr;
pval_t   = p.Results.pval;
verbose  = p.Results.verbose;
w12mode  = string(p.Results.w12mode);
alpha_   = p.Results.alpha;

X_s = double(full(X_s));
X_t = double(full(X_t));
g   = string(g(:));

ng = size(X_s, 1);
if size(X_t, 1) ~= ng || numel(g) ~= ng
    error('xctmain_nn: X_s, X_t and g must all have the same number of rows.');
end
% if ng > 3000
%     warning('xctmain_nn:largeGeneSet', ...
%         ['%d genes → L matrix will be %d × %d (%.0f MB single). ' ...
%          'Consider subsetting to top HVGs to reduce memory.'], ...
%         ng, 2*ng, 2*ng, (2*ng)^2*4/1e6);
% end

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
cfg = struct('n_dim', n_dim, 'mu', mu_, 'ncomp', ncomp, 'grn_q', grn_q, ...
    'n_steps', n_steps, 'lr', lr, 'pval', pval_t, 'verbose', verbose, ...
    'w12mode', w12mode, 'alpha', alpha_, 'grnoffset', p.Results.grnoffset, ...
    'useparallel', logical(p.Results.useparallel));

[T1, E1] = i_align_nn(X_s, X_t, g, lig_db, rec_db, cfg);

if twosided
    [T2, E2] = i_align_nn(X_t, X_s, g, lig_db, rec_db, cfg);
    T  = {T1, T2};
    emb = {E1, E2};
else
    T = T1;
    emb = E1;
end

end % xctmain_nn


%% ── CORE ALIGNMENT ───────────────────────────────────────────────────────

function [T, emb] = i_align_nn(X_s, X_t, g, lig_db, rec_db, cfg)
% I_ALIGN_NN  Neural-network manifold alignment for one direction.

n_dim   = cfg.n_dim;
mu_     = cfg.mu;
ncomp   = cfg.ncomp;
grn_q   = cfg.grn_q;
n_steps = cfg.n_steps;
lr      = cfg.lr;
pval_t  = cfg.pval;
verbose = cfg.verbose;
w12mode = cfg.w12mode;
alpha_  = cfg.alpha;
grnoffset = cfg.grnoffset;
useparallel = cfg.useparallel;

ng        = size(X_s, 1);
g_up      = upper(g);

% ── 1. Within-type PCNet GRNs ────────────────────────────────────────────
if verbose, fprintf('[xctmain_nn]   Building PCNet GRNs ...\n'); end
W11 = ten.i_xctgrn(X_s, ncomp, grn_q, verbose, useparallel);   % sparse ng × ng
W22 = ten.i_xctgrn(X_t, ncomp, grn_q, verbose, useparallel);   % sparse ng × ng

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
emb = struct('P_s', [], 'P_t', [], 'li_idx', li_idx, 'ri_idx', ri_idx);
if n_hit == 0
    warning('xctmain_nn:noPairs', 'No L-R pairs found. Returning empty table.');
    T = table(); return;
end
if verbose, fprintf('[xctmain_nn]   %d L-R pairs matched.\n', n_hit); end

% Correspondence block. "outer" reproduces core.py:_build_w with its default
% query_DB=None, in which the L-R database does not enter the alignment at all
% and is applied only when candidates are ranked below.
W12 = ten.i_xctw12(X_s, X_t, ng, mode=w12mode, alpha=alpha_, ...
    li_idx=li_idx, ri_idx=ri_idx, lognorm=false);

% ── 3. Block weight matrix W and graph Laplacian L ───────────────────────
%   The offset is added to EVERY GRN element, as core.py:415 does. An earlier
%   comment here claimed this file followed the reference "exactly - no +1
%   shift"; that was backwards, and the omission is what left the graph
%   disconnected. See TEN.I_XCTBLOCK.
W = ten.i_xctblock(W11, W22, W12, mu_, grnoffset);

% ── 4-6. Train the networks and extract the embedding ────────────────────
%   Delegated to ten.i_nnembed so that this function and ten.sctenifoldxct
%   share one implementation of the published training loop.
if verbose
    fprintf('[xctmain_nn]   Training (%d steps, lr=%.4f) ...\n', n_steps, lr);
end

[P_s, P_t] = ten.i_nnembed(W, ng, n_dim, X_s, X_t, ...
    n_steps=n_steps, lr=lr, verbose=verbose);

% Returned so ten.xct.xctmain2_nn can build its chi-square statistic over all
% gene pairs, which the candidate-only output table cannot support.
emb.P_s = P_s;
emb.P_t = P_t;

% ── 7. Candidate L-R distances ────────────────────────────────────────────
cand_d = zeros(n_hit, 1);
for k = 1:n_hit
    d_k      = P_s(li_idx(k),:) - P_t(ri_idx(k),:);
    cand_d(k) = sqrt(d_k * d_k');
end

% ── 8. Left-tail null test ───────────────────────────────────────────────
% Exhaustive over non-candidate pairs, as stat.py::null_test. Small p-value
% means the candidate distance is unusually small. See TEN.I_XCTNULL for why
% this replaced a sampled null.
p_vals = ten.i_xctnull(P_s, P_t, li_idx, ri_idx, W12);

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


%% ── SHARED HELPERS (duplicated from xctmain for self-containment) ────────

function i_rejectcorrthr(args)
% I_REJECTCORRTHR  Explain the removal of the Pearson GRN option.
%   The 'corr_thr' argument selected a cut-off for a Pearson co-expression
%   proxy that no longer exists, so silently ignoring it would leave callers
%   believing they had configured a network they had not.

isname = cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, 'corr_thr'), args);
if any(isname)
    error("TEN:XCT:XCTMAIN_NN:CorrThrRemoved", ...
        "'corr_thr' is no longer supported. Within-type networks are now " + ...
        "built with PCNet rather than thresholded Pearson co-expression, " + ...
        "to match the published Python implementation. Use 'ncomp' to set " + ...
        "the number of principal components, or 'grn_q' to filter edges.");
end

end % i_rejectcorrthr


function X = i_lognorm(X)
% I_LOGNORM  Library-size normalisation + log1p.
cs = sum(X, 1);
cs(cs == 0) = 1;
X  = log1p(X ./ cs .* median(cs));
end % i_lognorm
