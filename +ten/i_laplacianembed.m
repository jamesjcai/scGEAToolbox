function [P_s, P_t, n_used] = i_laplacianembed(W, ng, n_dim, verbose)
%I_LAPLACIANEMBED  Spectral embedding of a two-block manifold-alignment graph.
%   [P_S, P_T] = TEN.I_LAPLACIANEMBED(W, NG, N_DIM) builds the unnormalized
%   graph Laplacian of the 2*NG-by-2*NG block weight matrix W and returns the
%   source and target gene coordinates given by its N_DIM smallest non-trivial
%   eigenvectors.
%
%   INPUTS:
%     W       - 2*NG-by-2*NG symmetric block weight matrix
%               [W_source, mu*W12; mu*W12', W_target]
%     ng      - number of genes (half the side length of W)
%     n_dim   - requested embedding dimensions
%     verbose - print progress (default false)
%
%   OUTPUTS:
%     P_s    - NG-by-N_USED source gene coordinates
%     P_t    - NG-by-N_USED target gene coordinates
%     n_used - embedding dimensions actually used (<= N_DIM)
%
%   Two numerical details matter here, and getting either wrong yields an
%   embedding of all zeros rather than an error:
%
%   1. Convergence. EIGS with 'smallestreal' routinely fails to converge on a
%      Laplacian of this size and returns NaN eigenvalues. Shift-invert about a
%      small negative sigma converges reliably and is several times faster.
%      Sigma must be nonzero because a Laplacian is exactly singular, so
%      'smallestabs' cannot factorize it.
%
%   2. Null-space multiplicity. The multiplicity of eigenvalue 0 equals the
%      number of connected components of the graph, and a thresholded GRN is
%      routinely disconnected. Requesting only N_DIM+4 eigenvalues can return
%      nothing but null space, leaving no non-trivial eigenvectors at all. The
%      component count is therefore measured first and the request sized past
%      it.
%
%   If no non-trivial eigenvector can be found the function errors rather than
%   returning a degenerate all-zero embedding, which would otherwise propagate
%   as a full table of zero distances and p-values of 1.
%
% see also: TEN.SCTENIFOLDXCT, TEN.SCTENIFOLDXCT2, TEN.XCT.XCTMAIN

arguments
    W {mustBeNumeric}
    ng (1, 1) double
    n_dim (1, 1) double
    verbose (1, 1) logical = false
end

n = 2*ng;

% ── Unnormalized Laplacian (|W| degree, so signed edges still give d >= 0) ──
d = full(sum(abs(W), 2));
L = spdiags(d, 0, n, n) - W;

% ── Size the request past the null space ──────────────────────────────────
ncomp = i_ncomponents(W);
n_ev = min(n_dim + ncomp + 2, n - 2);
n_ev = max(n_ev, min(3, n - 2));
if verbose
    fprintf('[i_laplacianembed] %d components; requesting %d eigenvalues.\n', ...
        ncomp, n_ev);
end

% ── Shift-invert eigensolve, with a plain fallback ────────────────────────
sigma = -1e-4;
try
    [V, D_ev] = eigs(L, n_ev, sigma);
catch
    % Fall back to the direct formulation if the shifted factorization fails
    [V, D_ev] = eigs(L, n_ev, 'smallestreal', struct('tol', 1e-6));
end

ev = real(diag(D_ev));
V = real(V);
ok = ~isnan(ev);
ev = ev(ok);
V = V(:, ok);
if isempty(ev)
    error("TEN:I_LAPLACIANEMBED:NoConvergence", ...
        "The eigensolver returned no converged eigenvalues for a %d-by-%d " + ...
        "Laplacian. Reduce n_dim, or filter the expression matrix to fewer " + ...
        "genes so the gene network is better conditioned.", n, n);
end

[ev, ord] = sort(ev, 'ascend');
V = V(:, ord);

% ── Drop trivial (null-space) eigenvectors, relative to the spectrum size ──
tol0 = 1e-8*max(1, max(abs(ev)));
keep = ev >= tol0;
V = V(:, keep);

if isempty(V)
    error("TEN:I_LAPLACIANEMBED:NoNonTrivial", ...
        "All %d computed eigenvalues lie in the null space, so no embedding " + ...
        "exists. The alignment graph has %d connected components; lowering the " + ...
        "GRN edge threshold or filtering low-expression genes will connect it.", ...
        numel(ev), ncomp);
end

n_used = min(n_dim, size(V, 2));
if n_used < n_dim && verbose
    warning("TEN:I_LAPLACIANEMBED:DimReduced", ...
        "Only %d non-trivial eigenvectors available; n_dim reduced from %d.", ...
        n_used, n_dim);
end
V = V(:, 1:n_used);

P_s = V(1:ng, :);
P_t = V(ng+1:end, :);

end


%% ---- number of connected components of the graph underlying W ----
function ncomp = i_ncomponents(W)
A = abs(W) > 0;
A = A | A';
try
    ncomp = max(conncomp(graph(A, "omitselfloops")));
catch
    % graph/conncomp unavailable; assume connected and let the caller proceed
    ncomp = 1;
end
end
