function [T, grns] = i_xctcore(X, g, ctype, celltype1, celltype2, twosided, cfg)
%I_XCTCORE  Shared engine for the scTenifoldXct manifold-alignment entry points.
%   T = TEN.I_XCTCORE(X, G, CTYPE, CELLTYPE1, CELLTYPE2, TWOSIDED, CFG) runs the
%   full alignment: two GRNs, a cross-type correspondence block, a spectral or
%   neural embedding, and a left-tail null test on the candidate pair distances.
%
%   This exists so that TEN.SCTENIFOLDXCT and TEN.SCTENIFOLDXCT_GLYCO cannot
%   drift apart. Everything that is common to both lives here; the only thing
%   either may vary is the content of the correspondence block, supplied through
%   cfg.corrfcn. This function has no knowledge of glycobiology.
%
%   INPUTS:
%     X         - genes-by-cells expression matrix, already subset to the two
%                 cell types and stored as the caller holds it (no transform is
%                 applied before the GRN step)
%     g         - gene names, one per row of X
%     ctype     - cell-type label per column of X
%     CELLTYPE1 - source cell type label
%     CELLTYPE2 - target cell type label
%     TWOSIDED  - logical; when true both directions are run and T is {T1, T2}
%
%   CFG FIELDS:
%     n_dim      - embedding dimensions
%     mu         - cross-cell-type correspondence weight
%     pval       - p-value threshold applied to the output table
%     verbose    - print progress
%     w12mode    - "lr" (cognate pairs) or "outer" (reference outer product)
%     alpha      - mean/variance blend for w12mode="outer"
%     grnoffset  - constant added to every GRN element before assembly; see
%                  TEN.I_XCTBLOCK
%     slv        - struct with .name ("spectral"|"nn"), .n_steps, .lr, .seed
%     tag        - string used to prefix progress messages
%     corrfcn    - [] for the published weight-1 correspondence, or a handle
%                     [li, ri, w, module, channel] = corrfcn(g, li, ri, src, tgt)
%                  which may reweight the matched pairs and append new ones.
%                  Called once per direction, with src/tgt naming that direction.
%     provenance - logical; when true the output gains channel, glyco_module and
%                  glyco_weight columns describing the correspondence actually
%                  used. Left false the schema is exactly the published one.
%     grn1, grn2 - [] to build each cell type's GRN via TEN.I_XCTGRN as usual,
%                  or an ng-by-ng precomputed adjacency (same gene order as g)
%                  to use for CELLTYPE1/CELLTYPE2 instead, skipping the
%                  net.pcrnet build for that side. See TEN.I_XCTGRN's
%                  'precomputed' option for what still happens to it (scale,
%                  optional filter, symmetrize) and why substituting a network
%                  built by a different method is a real methodological choice.
%     grn1_processed, grn2_processed - logical; true if the matching grn1/
%                  grn2 already went through TEN.I_XCTGRN's scale/filter/
%                  symmetrize (e.g. it is a grns.A_s/A_t this function itself
%                  returned earlier) and that should not be done again - see
%                  TEN.I_XCTGRN's 'processed' option for why a second pass is
%                  not simply a no-op. Ignored when the matching grn1/grn2 is
%                  empty. Default false.
%     candidates - "database" (default) restricts the output to L-R database
%                  matches, as published. "all" instead ranks EVERY gene pair
%                  by embedding distance and returns the closest cfg.topN,
%                  regardless of whether either gene is a known ligand or
%                  receptor - a discovery mode for finding co-regulated pairs
%                  the database does not already know about. The output
%                  schema differs: no p_value/null test (there is no separate
%                  background left once every pair is a candidate), and a
%                  percentile column instead, plus is_known_lr marking which
%                  of the top pairs happen to already be in the database. A
%                  fully database-free run also needs w12mode="outer", since
%                  w12mode="lr" still uses the database to build the
%                  alignment itself even when candidates="all" only changes
%                  what gets reported afterward.
%
%   OUTPUT:
%     T    - table with columns ligand, receptor, dist, correspondence,
%            p_value, plus the three provenance columns when cfg.provenance
%            is true. When TWOSIDED is true, T is {T1, T2} for the two
%            directions.
%     grns - requested via a second output argument (nargout > 1); struct
%            with fields A_s, A_t (the ng-by-ng adjacency actually used for
%            CELLTYPE1/CELLTYPE2 - built fresh or cfg.grn1/grn2 passed
%            through TEN.I_XCTGRN's scale/filter/symmetrize either way) and
%            genes (g, for traceability). One pair regardless of TWOSIDED:
%            both directions align the same two networks, just swapped.
%            Building this struct costs nothing extra - A_s/A_t already exist
%            by the time this could be returned - so it is always populated
%            when asked for, never gated behind its own cfg flag.
%
% see also: TEN.SCTENIFOLDXCT, TEN.SCTENIFOLDXCT_GLYCO, TEN.I_ALIGNEMBED,
%           TEN.I_XCTW12

tag = cfg.tag;
verbose = cfg.verbose;

is_s = ctype == celltype1;
is_t = ctype == celltype2;

X_s = X(:, is_s);   % genes x source cells
X_t = X(:, is_t);   % genes x target cells

if verbose
    fprintf('[%s] %s: %d genes x %d cells\n', tag, celltype1, ...
        size(X_s, 1), size(X_s, 2));
    fprintf('[%s] %s: %d genes x %d cells\n', tag, celltype2, ...
        size(X_t, 1), size(X_t, 2));
end

% -- Ligand-receptor database --------------------------------------------
[lig_db, rec_db] = i_loadlrdb();
if verbose
    fprintf('[%s] L-R database: %d pairs loaded.\n', tag, numel(lig_db));
end

% -- GRNs -----------------------------------------------------------------
% GRNs via the shared builder. This function's settings differ from the
% xctmain path's - ncomp=3 with a 0.75 edge filter, against the reference's
% ncomp=5 unfiltered - but the code is now the same, so the two cannot drift
% apart in anything except their arguments. TEN.I_XCTGRN owns the pcrnet call,
% the max-normalisation, the filtering and the symmetrisation, and records why
% parallel and fastersvd are both left off.
if verbose
    if isempty(cfg.grn1)
        fprintf('[%s] Building GRN: %s ...\n', tag, celltype1);
    else
        fprintf('[%s] Using precomputed GRN: %s (skipping pcrnet build)\n', tag, celltype1);
    end
end
A_s = ten.i_xctgrn(X_s, 3, 0.75, false, cfg.useparallel, ...
    precomputed=cfg.grn1, processed=cfg.grn1_processed);

if verbose
    if isempty(cfg.grn2)
        fprintf('[%s] Building GRN: %s ...\n', tag, celltype2);
    else
        fprintf('[%s] Using precomputed GRN: %s (skipping pcrnet build)\n', tag, celltype2);
    end
end
A_t = ten.i_xctgrn(X_t, 3, 0.75, false, cfg.useparallel, ...
    precomputed=cfg.grn2, processed=cfg.grn2_processed);

% -- Alignment, one direction at a time -----------------------------------
if verbose
    fprintf('[%s] Aligning %s -> %s ...\n', tag, celltype1, celltype2);
end
T1 = i_xct(X_s, X_t, g, A_s, A_t, lig_db, rec_db, ...
    string(celltype1), string(celltype2), cfg);

if twosided
    if verbose
        fprintf('[%s] Aligning %s -> %s ...\n', tag, celltype2, celltype1);
    end
    T2 = i_xct(X_t, X_s, g, A_t, A_s, lig_db, rec_db, ...
        string(celltype2), string(celltype1), cfg);
    T = {T1, T2};
else
    T = T1;
end

if nargout > 1
    grns = struct('A_s', A_s, 'A_t', A_t, 'genes', string(g(:)), ...
        'celltype1', string(celltype1), 'celltype2', string(celltype2));
end

end % i_xctcore


%% ---- manifold alignment for one direction (source -> target) ----
function T = i_xct(X_s, X_t, g, A_s, A_t, lig_db, rec_db, src, tgt, cfg)

n_dim = cfg.n_dim;
pval_t = cfg.pval;
verbose = cfg.verbose;
tag = cfg.tag;

ng = size(X_s, 1);          % number of genes (same for source and target)

% -- Match database pairs against the gene list ---------------------------
% Matched unconditionally even in candidates="all": w12mode="lr" needs it
% for the correspondence block regardless of candidate scope, and even a
% w12mode="outer" discovery run uses it to flag which top pairs happen to
% already be known (is_known_lr), so the two are not fully independent
% pieces of information even when the OUTPUT does not stop at the database.
[li_idx, ri_idx] = i_matchpairs(g, lig_db, rec_db);
n_valid = numel(li_idx);

if cfg.candidates == "database" && n_valid == 0
    warning('sctenifoldxct:noPairs', ...
        'No L-R pairs found in gene list. Returning empty table.');
    T = table();
    return;
end
if cfg.w12mode == "lr" && n_valid == 0
    warning('sctenifoldxct:noPairs', ...
        ['No L-R pairs found in gene list; w12mode="lr" needs at least one pair ' ...
        'to build the correspondence block. Returning empty table.']);
    T = table();
    return;
end
if verbose && n_valid > 0
    fprintf('[%s]   %d L-R pairs matched in data.\n', tag, n_valid);
end

% -- Correspondence weights -----------------------------------------------
% Default: every pair gets weight 1, exactly as published. A corrfcn may
% reweight these and append further correspondences.
lr_w = ones(n_valid, 1);
lr_module = repmat("", n_valid, 1);
channel = repmat("protein_LR", n_valid, 1);

if ~isempty(cfg.corrfcn)
    [li_idx, ri_idx, lr_w, lr_module, channel] = ...
        cfg.corrfcn(g, li_idx, ri_idx, src, tgt);
end

% -- Correspondence block -------------------------------------------------
% "lr" aggregates duplicate cognate pairs by summing their weights. "outer"
% instead reproduces the reference correspondence, a dense expression outer
% product over all gene pairs in which the L-R database plays no part; the
% matched pairs are then used only to name candidates and measure distances.
% See TEN.I_XCTW12.
if cfg.w12mode == "outer"
    W12 = ten.i_xctw12(X_s, X_t, ng, mode="outer", alpha=cfg.alpha);
else
    W12 = sparse(li_idx, ri_idx, lr_w, ng, ng);
end

% -- Block weight matrix --------------------------------------------------
% Assembly, the GRN offset and the mu rescaling all live in TEN.I_XCTBLOCK.
% Because of that renormalisation the total cross-type coupling is invariant to
% the correspondence weights; only the relative pull between pairs changes, so
% mu does not need retuning when a corrfcn reweights them.
W = ten.i_xctblock(0.5*(A_s + A_s'), 0.5*(A_t + A_t'), W12, ...
    cfg.mu, cfg.grnoffset);

% -- Embedding ------------------------------------------------------------
% Delegated so both methods share one interface and the eigensolver's
% convergence and null-space handling stay in one place. See TEN.I_ALIGNEMBED
% for why the two are different methods rather than interchangeable solvers.
if cfg.slv.name == "nn"
    % The network consumes expression profiles directly through sigmoid
    % layers, so its inputs are log-normalised here to keep them off the
    % saturated tails. The spectral path is unaffected: it reads only W.
    [P_s, P_t] = ten.i_alignembed(W, ng, n_dim, solver="nn", ...
        X_s=i_lognorm(X_s), X_t=i_lognorm(X_t), ...
        n_steps=cfg.slv.n_steps, lr=cfg.slv.lr, seed=cfg.slv.seed, ...
        verbose=verbose);
else
    [P_s, P_t] = ten.i_alignembed(W, ng, n_dim, solver="spectral", ...
        verbose=verbose);
end

if cfg.candidates == "all"
    T = i_xctdiscover(P_s, P_t, g, li_idx, ri_idx, cfg.topN, tag, verbose);
    return
end

% -- Candidate distances --------------------------------------------------
n_cand = numel(li_idx);
cand_d = zeros(n_cand, 1, 'single');
for k = 1:n_cand
    diff_k = P_s(li_idx(k), :) - P_t(ri_idx(k), :);
    cand_d(k) = sqrt(sum(diff_k.^2));
end

% -- Left-tail null test --------------------------------------------------
% Delegated to TEN.I_XCTNULL, which scores against the exhaustive set of
% non-candidate pairs exactly as stat.py::null_test does. This used to draw a
% random sample, which made the p-values RNG-dependent; see that function's
% help. Small p_value -> candidate distance unusually small -> strong
% interaction.
p_vals = ten.i_xctnull(P_s, P_t, li_idx, ri_idx, W12);

% -- Output table ---------------------------------------------------------
lig_out = g(li_idx);
rec_out = g(ri_idx);

% Per-pair W12 weight actually used, summed over duplicate entries if any.
w12_vals = full(W12(sub2ind([ng, ng], li_idx, ri_idx)));

T = table(lig_out, rec_out, double(cand_d), w12_vals, p_vals, ...
    'VariableNames', {'ligand', 'receptor', 'dist', 'correspondence', 'p_value'});

% Provenance columns are appended only on request, so the default output
% schema is exactly the published one.
if cfg.provenance
    T.channel = channel;
    T.glyco_module = lr_module;
    T.glyco_weight = lr_w;
end

T = T(T.p_value <= pval_t, :);
T = sortrows(T, 'dist', 'ascend');

if verbose
    fprintf('[%s]   %d significant L-R pairs (pval <= %.2f).\n', ...
        tag, height(T), pval_t);
end

end % i_xct


%% ---- discovery mode: rank every gene pair by embedding distance ----
function T = i_xctdiscover(P_s, P_t, g, li_idx, ri_idx, topN, tag, verbose)
%I_XCTDISCOVER  Every gene in the source embedded against every gene in the
%target, ranked by distance, with no L-R database restriction on the
%output. There is no longer a separate "non-candidate" population to test
%against (every pair IS a candidate), so instead of TEN.I_XCTNULL's p-value
%this reports each pair's percentile among all pairs directly - equivalent
%math (a small value still means "closer than almost everything else"), just
%computed from the same rank the sort already produced rather than a
%permutation test against a background that no longer exists.
ng = size(P_s, 1);
D = pdist2(P_s, P_t);   % ng-by-ng; D(i,j) = dist(gene i source, gene j target)
D(1:ng+1:end) = Inf;    % a gene is not its own ligand/receptor; excluded, not just deprioritized

[sortedD, lin] = sort(D(:), 'ascend');
nPairs = ng * (ng - 1);   % off-diagonal pairs only
keep = 1:min(topN, nPairs);
lin = lin(keep);
[si, ti] = ind2sub([ng, ng], lin);

known = false(numel(lin), 1);
if ~isempty(li_idx)
    known = ismember(lin, sub2ind([ng, ng], li_idx, ri_idx));
end

T = table(g(si), g(ti), double(sortedD(keep)), keep(:)/nPairs, known, ...
    'VariableNames', {'ligand', 'receptor', 'dist', 'percentile', 'is_known_lr'});

if verbose
    fprintf('[%s]   candidates=all: top %d of %d gene pairs by distance (%d already in the L-R database).\n', ...
        tag, height(T), nPairs, sum(T.is_known_lr));
end
end % i_xctdiscover


%% ---- indices of database pairs present in the gene list ----
function [li_idx, ri_idx] = i_matchpairs(g, lig_db, rec_db)
g_up = upper(string(g(:)));
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
end % i_matchpairs


%% ---- built-in ligand-receptor database ----
function [lig_db, rec_db] = i_loadlrdb()
pw1 = fileparts(mfilename('fullpath'));
lr_mat = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.mat');
lr_txt = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.txt');

if exist(lr_mat, 'file')
    db = load(lr_mat, 'ligand', 'receptor');
    lig_db = upper(string(db.ligand(:)));
    rec_db = upper(string(db.receptor(:)));
else
    % Fallback: read from text file (tab-delimited, cols 3 and 5)
    T_lr = readtable(lr_txt, 'FileType', 'text', 'Delimiter', '\t');
    lig_db = upper(string(T_lr{:, 3}));
    rec_db = upper(string(T_lr{:, 5}));
end
end % i_loadlrdb


%% ---- library-size normalisation followed by log1p ----
function X = i_lognorm(X)
% Applied only to the neural solver's network inputs; the GRN step keeps the
% caller's original matrix.

X = double(X);
col_sums = sum(X, 1);
col_sums(col_sums == 0) = 1;    % avoid /0 for empty cells
X = X./col_sums.*median(col_sums);
X = log1p(X);
end % i_lognorm
