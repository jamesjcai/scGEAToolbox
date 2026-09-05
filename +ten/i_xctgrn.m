function A = i_xctgrn(X, ncomp, q, verbose, useparallel, opts)
%I_XCTGRN  Within-type GRN for the scTenifoldXct entry points.
%   A = TEN.I_XCTGRN(X) builds a principal-component regression network from
%   the genes-by-cells matrix X and returns it as a sparse, symmetric,
%   max-normalized adjacency matrix.
%
%   A = TEN.I_XCTGRN(X, NCOMP, Q, VERBOSE) sets the number of principal
%   components, an optional edge-filtering quantile, and progress printing.
%
%   The defaults reproduce the configuration the published Python
%   scTenifoldXct uses. core.py calls make_pcNet(nComp=5), and pcNet.py then
%   applies scale=True, q=0. and symmetric=True, so the network is scaled by
%   its largest absolute entry, left unthresholded, and symmetrized last. The
%   ordering matters: symmetrizing before scaling would change the maximum.
%
%   Two differences from pcNet.py remain and are not worth reimplementing:
%   net.pcrnet z-scores genes, whereas pcNet.py L2-normalizes columns without
%   centering; and Q here delegates to TEN.E_FILTADJC, which takes a quantile
%   of the nonzero entries after rescaling by the top-10 mean, rather than a
%   percentile of all entries. Q defaults to 0 (no filtering), so the second
%   difference is inert unless a caller asks for it.
%
%   INPUTS:
%     X       - genes-by-cells expression matrix, log-normalized
%     ncomp   - principal components (default 5, as in core.py)
%     q       - edge-filtering quantile in [0, 1); 0 disables (default 0)
%     verbose - print timing (default false)
%
%   NAME-VALUE ARGUMENTS:
%     symmetrize - apply 0.5*(A+A') last, as pcNet.py's symmetric=True does
%                  (default true). FALSE for scTenifoldNet and scTenifoldKnk,
%                  which need the network directed; see below.
%     scale      - normalize by max(abs(A)) (default true). FALSE for
%                  TEN.I_NCT, which filters and saves the raw network itself.
%     fastersvd  - swap svds for lmsvd inside net.pcrnet (default false).
%                  Slower, and non-deterministic; see below.
%     precomputed - an ng-by-ng adjacency to use in place of net.pcrnet's
%                  output, e.g. a network built by another tool for the same
%                  genes and cell population. Its diagonal is zeroed (a
%                  precomputed network may carry self-terms net.pcrnet's
%                  regression never produces) and it is then scaled/filtered/
%                  symmetrized exactly as a freshly-built one would be, so
%                  cfg.mu and cfg.grnoffset stay meaningful downstream. X is
%                  still required, for its size only - X is genes-by-cells,
%                  precomputed must be genes-by-genes for the same genes.
%                  Substituting a network built by a different method (a
%                  different ncomp, a bootstrapped/denoised one, ...) changes
%                  what is being tested, not just how it was computed; that is
%                  a deliberate choice for the caller to make, not a hidden one.
%     processed  - the PRECOMPUTED network already went through this exact
%                  scale/filter/symmetrize pipeline (e.g. it is the A_s/A_t
%                  TEN.SCTENIFOLDXCT's grns output returned from an earlier
%                  call) - skip doing it again and use PRECOMPUTED as-is
%                  (default false). This exists because the steps are not all
%                  idempotent: Q-filtering an already-filtered network filters
%                  it again relative to its now-smaller nonzero set, so a
%                  second pass at the same Q removes MORE edges, not the same
%                  ones - net.pcrnet(X,3) -> filter once -> filter again is a
%                  different, sparser network than net.pcrnet(X,3) -> filter
%                  once, not the same result recomputed. Scale and symmetrize
%                  are themselves idempotent, so this only matters because Q
%                  is not; it is simplest to skip all three together whenever
%                  the input has already been through them.
%
%   OUTPUT:
%     A - ng-by-ng adjacency matrix; sparse and symmetric by default
%
%   SHARED BY ALL THREE TOOLS. scTenifoldNet, scTenifoldKnk and scTenifoldXct
%   all build their within-type networks with PCNet, and all of them route
%   through this function, so the pcrnet invocation exists in one place. What
%   they do NOT share is the post-processing, and the differences are real
%   rather than accidental:
%
%     ten.xct.xctmain, xctmain_nn   ncomp=5, q=0,    scaled, symmetric
%     ten.sctenifoldxct, ...xct2    ncomp=3, q=0.75, scaled, symmetric
%     ten.i_nc  (scTenifoldNet/Knk) caller ncomp,    scaled, DIRECTED
%     ten.i_nct (scTenifoldNet)     caller ncomp,    raw,    DIRECTED
%
%   The directed cases are load-bearing: sctenifoldnet.m runs its tensor
%   decomposition on the raw directed stack and symmetrizes only afterwards,
%   immediately before manifold alignment. Symmetrizing here would move that
%   step earlier and change the algorithm.
%
%   Requires Statistics and Machine Learning Toolbox, via net.pcrnet's zscore.
%
% see also: NET.PCRNET, TEN.E_FILTADJC, TEN.I_NC, TEN.I_NCT,
%           TEN.XCT.XCTMAIN, TEN.XCT.XCTMAIN_NN

arguments
    X double
    ncomp (1, 1) double {mustBePositive} = 5
    q (1, 1) double {mustBeGreaterThanOrEqual(q, 0), mustBeLessThan(q, 1)} = 0
    verbose (1, 1) logical = false
    useparallel (1, 1) logical = false
    opts.symmetrize (1, 1) logical = true
    opts.fastersvd (1, 1) logical = false
    opts.scale (1, 1) logical = true
    opts.precomputed double = []
    opts.processed (1, 1) logical = false
end

ng = size(X, 1);
if ~isempty(opts.precomputed)
    if ~isequal(size(opts.precomputed), [ng, ng])
        error("TEN:I_XCTGRN:BadPrecomputedSize", ...
            "precomputed must be %d-by-%d (X has %d genes). Received %s.", ...
            ng, ng, ng, mat2str(size(opts.precomputed)));
    end
elseif ncomp < 2 || ncomp >= ng
    error("TEN:I_XCTGRN:BadNComp", ...
        "NCOMP must be at least 2 and fewer than the %d genes in X. " + ...
        "Received %g.", ng, ncomp);
end

% net.pcrnet regresses each gene on the principal components of the others, a
% loop over ng genes that it runs as a parfor when USEPARALLEL is set. A GPU
% takes precedence: pcrnet disables the parfor itself, because parfor and
% gpuArray cannot be combined.
%
% USEPARALLEL defaults to false because it is usually the slower choice, which
% is counterintuitive enough to be worth recording. The serial loop is ALREADY
% parallel: multithreaded BLAS spreads each SVD across every core, measured at
% 3.0x on 600 genes and 5.7x on 1856 (20 cores, 784 cells). A parfor fragments
% that same work and ends up fighting it - against serial with no pool, roughly
% break-even at 600 genes, 0.5x at 1200 and 0.3x at 1856.
%
% A Threads pool does not rescue it. Its workers report the full 20
% computational threads each, so 20 workers oversubscribe 20 cores 20-fold;
% measured 0.41x at 1856, no better than the process pool's 0.62x. Process
% workers have the opposite problem, one thread apiece. Either way the BLAS
% threading the serial path gets for free is the thing to preserve.
%
% pcrnet's FASTERSVD flag is likewise left off, but for reproducibility rather
% than for speed. It swaps svds for lmsvd, a block-Krylov solver.
%
% ON SPEED, THE HONEST ANSWER IS "ABOUT THE SAME". Early single-shot tic/toc
% runs suggested lmsvd was 1.7-2.5x slower, and that was wrong - this machine's
% wall clock varies by over 2x between runs on identical work (svds on 1200
% genes measured 6.8 s and 15.0 s in two sittings). Re-measured with timeit,
% per-call times run 0.53x to 1.13x with no consistent direction, and the whole
% pcrnet loop at 1200 genes comes out at 0.97x. Treat the two as comparable and
% do not quote a speedup either way. lmsvd converges in ~10 iterations against
% its maxit of 150, so it is not mistuned; it simply has no edge over ARPACK at
% the small r used here.
%
% THE REAL REASONS to prefer svds:
%   1. lmsvd is NON-DETERMINISTIC. Its default initial guess is randn(n,r), so
%      two identical calls differ by ~1e-9 where svds is bitwise identical.
%      That alone prevents scTenifoldNet and scTenifoldKnk from reproducing a
%      run exactly.
%   2. Accuracy degrades with r: relative difference 5e-9 at ncom=3 rising to
%      1.3e-5 at ncom=40.
%
% Both solvers compute the same subspace - |Vs'*Vl| has singular values
% [1 1 1] - so switching costs nothing in the result: the L-R ranking is
% unchanged, rank correlation +1.0000 with zero rank shift.
%
% fastersvd was turned off at every call site in the toolbox on 2026-07-20,
% including the scTenifoldNet and scTenifoldKnk paths that had used it since
% before this function existed.
if opts.processed && isempty(opts.precomputed)
    error("TEN:I_XCTGRN:ProcessedNeedsPrecomputed", ...
        "processed=true has nothing to skip processing ON without precomputed " + ...
        "also being given - a freshly built network still needs scale/filter/" + ...
        "symmetrize. Pass precomputed too, or leave processed at its default.");
end

t0 = tic;
if ~isempty(opts.precomputed)
    A = double(opts.precomputed);
    A(1:ng+1:end) = 0;   % net.pcrnet never produces self-terms; match that
else
    useGPU = pkg.i_usegpu(X);
    A = net.pcrnet(X, ncomp, opts.fastersvd, true, useparallel && ~useGPU, ...
        false, useGPU);
end

if ~opts.processed
    % scale=True: normalize by the largest absolute entry. Off for ten.i_nct,
    % which hands the raw network to its own filtering and saving step.
    if opts.scale
        mx = max(abs(A), [], "all");
        if mx > 0
            A = A./mx;
        end
    end

    % q > 0: optional edge filtering, off by default to match pcNet.py.
    if q > 0
        A = ten.e_filtadjc(A, q, false);
    end

    % symmetric=True, applied last as in pcNet.py.
    %
    % SYMMETRIZE=FALSE exists for scTenifoldNet and scTenifoldKnk, which must
    % keep the network directed: sctenifoldnet.m runs the tensor decomposition
    % on the raw asymmetric stack and only symmetrizes afterwards, right
    % before manifold alignment. Symmetrizing here would move that step
    % earlier and change the algorithm, so it is off for those callers and on
    % for every xct one.
    if opts.symmetrize
        A = sparse(0.5*(A + A.'));
    end
else
    A = sparse(A);   % match the shape (sparse) the processed branch produces
end

if verbose
    if isempty(opts.precomputed)
        fprintf('[i_xctgrn] %d genes, ncomp=%d, %.1f%% nonzero, %.1f s\n', ...
            ng, ncomp, 100*nnz(A)/ng^2, toc(t0));
    else
        fprintf('[i_xctgrn] %d genes, precomputed network, %.1f%% nonzero, %.1f s\n', ...
            ng, 100*nnz(A)/ng^2, toc(t0));
    end
end

end
