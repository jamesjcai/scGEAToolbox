function [W12, info] = i_xctw12(X_s, X_t, ng, opts)
%I_XCTW12  Cross-cell-type correspondence block for manifold alignment.
%   W12 = TEN.I_XCTW12(X_S, X_T, NG) builds the off-diagonal block of the
%   alignment graph. Two modes are available and they encode different ideas of
%   what a correspondence is.
%
%   mode="outer" (matches Python scTenifoldXct 0.2.0)
%     W12 = metric_s * metric_t', a dense outer product over EVERY gene pair,
%     with metric = (1-alpha)*mean^2 + alpha*var computed per gene on
%     log-normalised expression. The ligand-receptor database plays no part: the
%     correspondence says "both genes are expressed and variable", and the
%     database is applied downstream when candidate pairs are ranked. This
%     reproduces `core.py:_build_w` with its default `query_DB=None`, whose own
%     docstring reads "default None will not query the DB and use all pair-wise
%     corresponding scores".
%
%   mode="lr" (the toolbox's historical behaviour)
%     W12 is sparse with one entry per cognate ligand-receptor pair from the
%     database, weight 1 by default. The correspondence says "these two are a
%     known interacting pair".
%
%   These are not interchangeable. "outer" solves an expression-driven
%   alignment; "lr" solves a database-driven one. Choose by whether the goal is
%   parity with the published implementation or an L-R-structured analysis.
%
%   INPUTS:
%     X_s - NG-by-nCellsSource expression matrix (genes by cells, as stored)
%     X_t - NG-by-nCellsTarget expression matrix
%     ng  - number of genes
%
%   NAME-VALUE ARGUMENTS:
%     mode    - "outer" (default) or "lr"
%     alpha   - mean/variance blend for "outer", in [0,1] (default 0.5, as Python)
%     li_idx  - ligand row indices, required for "lr"
%     ri_idx  - receptor column indices, required for "lr"
%     weights - per-pair weights for "lr" (default all ones)
%     lognorm - log-normalise before computing metrics, "outer" only
%               (default true; set false if X is already log data)
%
%   OUTPUTS:
%     W12  - NG-by-NG correspondence block. Dense for "outer", sparse for "lr".
%     info - struct with .mode, .nnz, .dense, and for "outer" the metric vectors
%            .metric_s and .metric_t
%
%   MEMORY. "outer" is dense by construction, NG-by-NG, exactly as the reference
%   implementation is. At 2000 genes that is ~32 MB; at 20000 genes it is ~3.2 GB
%   and the assembled block matrix is four times that. Python has the same
%   ceiling and relies on the caller subsetting genes; so does this.
%
%   Variance uses the population convention (normalise by N) to match numpy's
%   `np.var` default of ddof=0. MATLAB's `var` default is the sample convention
%   and would not match.
%
% see also: TEN.SCTENIFOLDXCT, TEN.XCT.XCTMAIN, TEN.I_ALIGNEMBED

arguments
    X_s {mustBeNumeric}
    X_t {mustBeNumeric}
    ng (1, 1) double
    opts.mode (1, 1) string {mustBeMember(opts.mode, ["outer", "lr"])} = "outer"
    opts.alpha (1, 1) double {mustBeBetween(opts.alpha, 0, 1)} = 0.5
    opts.li_idx (:, 1) double = []
    opts.ri_idx (:, 1) double = []
    opts.weights (:, 1) double = []
    opts.lognorm (1, 1) logical = true
end

if opts.mode == "lr"
    if isempty(opts.li_idx) || isempty(opts.ri_idx)
        error("TEN:I_XCTW12:MissingPairs", ...
            "mode=""lr"" requires li_idx and ri_idx naming the cognate pairs.");
    end
    if numel(opts.li_idx) ~= numel(opts.ri_idx)
        error("TEN:I_XCTW12:PairSizeMismatch", ...
            "li_idx (%d) and ri_idx (%d) must have the same length.", ...
            numel(opts.li_idx), numel(opts.ri_idx));
    end
    w = opts.weights;
    if isempty(w)
        w = ones(numel(opts.li_idx), 1);
    elseif numel(w) ~= numel(opts.li_idx)
        error("TEN:I_XCTW12:WeightSizeMismatch", ...
            "weights (%d) must have one entry per pair (%d).", ...
            numel(w), numel(opts.li_idx));
    end

    W12 = sparse(opts.li_idx, opts.ri_idx, w, ng, ng);
    info = struct("mode", "lr", "nnz", nnz(W12), "dense", false);
    return;
end

% ── "outer": per-gene mean/variance metrics, then their outer product ─────
if size(X_s, 1) ~= ng || size(X_t, 1) ~= ng
    error("TEN:I_XCTW12:BadSizes", ...
        "X_S and X_T must both have %d rows (one per gene).", ng);
end

Xs = double(full(X_s));
Xt = double(full(X_t));
if opts.lognorm
    % core.py:_get_metric refuses integer input outright ("require log data"),
    % so the reference always sees log-normalised expression here.
    Xs = i_lognorm(Xs);
    Xt = i_lognorm(Xt);
end

metric_s = i_metric(Xs, opts.alpha);
metric_t = i_metric(Xt, opts.alpha);

W12 = metric_s*metric_t';       % dense ng x ng, as in core.py:414

info = struct("mode", "outer", "nnz", nnz(W12), "dense", true, ...
    "metric_s", metric_s, "metric_t", metric_t);

end


%% ---- (1-alpha)*mean^2 + alpha*var, per gene, population variance ----
function m = i_metric(X, alpha)
mu = mean(X, 2);
v = var(X, 1, 2);      % normalise by N, matching numpy's ddof=0 default
m = (1 - alpha).*mu.^2 + alpha.*v;
end


%% ---- library-size normalisation followed by log1p ----
function X = i_lognorm(X)
cs = sum(X, 1);
cs(cs == 0) = 1;
X = log1p(X./cs.*median(cs));
end
