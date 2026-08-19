function [P_s, P_t, info] = i_alignembed(W, ng, n_dim, opts)
%I_ALIGNEMBED  Embed a two-block manifold-alignment graph, by either method.
%   [P_S, P_T] = TEN.I_ALIGNEMBED(W, NG, N_DIM) embeds the genes of a
%   2*NG-by-2*NG block alignment graph using the spectral method.
%
%   [P_S, P_T] = TEN.I_ALIGNEMBED(..., solver="nn", X_s=XS, X_t=XT) instead
%   trains the neural-network alignment of the published scTenifoldXct.
%
%   Both minimise trace(P'*L*P), but they are different methods rather than two
%   solvers for one problem, because they optimise over different parameter
%   spaces:
%
%     "spectral" - P is free; the constrained minimum is the bottom
%                  eigenvectors of L, obtained in closed form. Global optimum,
%                  no toolbox dependency, deterministic.
%     "nn"       - P is the retraction of two MLPs applied to the expression
%                  matrices; the network weights are trained. The attained loss
%                  is generally higher, because the embedding is confined to
%                  the image of that function class. That restriction is a
%                  prior tying embedding to expression similarity.
%
%   Use "spectral" for any analytical claim. The "nn" route has two legitimate
%   uses - reproducing the published implementation, and benchmarking against
%   it - but it is NOT a robustness check on a prior encoded in W. That the
%   prior enters through W, upstream of both, does not make agreement across
%   solvers informative: it was measured and did not hold. On the glyco prior
%   the spectral route separated modulated from unmodulated pairs by a factor
%   of 38.9 while the neural route gave 0.64, and per-pair displacements
%   correlated at -0.055 between the two. A small change to W is absorbed into
%   a globally different but comparably-good network solution, and the
%   resulting run-to-run variation swamps the prior's signal. See
%   pipelines/glyco_xct/Notes.md section 6b.
%
%   INPUTS:
%     W     - 2*NG-by-2*NG symmetric block weight matrix
%     ng    - number of genes (half the side length of W)
%     n_dim - embedding dimensions
%
%   NAME-VALUE ARGUMENTS:
%     solver  - "spectral" (default) or "nn"
%     X_s     - NG-by-nCellsSource expression matrix; required for "nn"
%     X_t     - NG-by-nCellsTarget expression matrix; required for "nn"
%     n_steps  - Adam iterations, "nn" only (default 1000)
%     lr       - Adam learning rate, "nn" only (default 0.01)
%     seed     - RNG seed, "nn" only (default 0)
%     gradient - "riemannian" (default) or "euclidean", "nn" only. See
%                TEN.I_NNEMBED; "euclidean" reproduces pre-fix runs and is
%                measurably worse on the shared objective.
%     degree   - "signed" (default, as the reference) or "abs", "nn" only.
%                See the DEGREE CONVENTION section of TEN.I_NNEMBED.
%     verbose - print progress (default false)
%
%   OUTPUTS:
%     P_s, P_t - NG-by-K source and target coordinates
%     info     - struct; .solver always present, plus the diagnostics of
%                whichever method ran (see TEN.I_NNEMBED for the "nn" fields)
%
% see also: TEN.I_LAPLACIANEMBED, TEN.I_NNEMBED, TEN.SCTENIFOLDXCT

arguments
    W {mustBeNumeric}
    ng (1, 1) double
    n_dim (1, 1) double
    opts.solver (1, 1) string {mustBeMember(opts.solver, ["spectral", "nn"])} = "spectral"
    opts.X_s {mustBeNumeric} = []
    opts.X_t {mustBeNumeric} = []
    opts.n_steps (1, 1) double = 1000
    opts.lr (1, 1) double = 0.01
    opts.seed (1, 1) double = 0
    opts.gradient (1, 1) string ...
        {mustBeMember(opts.gradient, ["riemannian", "euclidean"])} = "riemannian"
    opts.degree (1, 1) string ...
        {mustBeMember(opts.degree, ["signed", "abs"])} = "signed"
    opts.verbose (1, 1) logical = false
end

switch opts.solver
    case "spectral"
        [P_s, P_t, n_used] = ten.i_laplacianembed(W, ng, n_dim, opts.verbose);
        info = struct("solver", "spectral", "n_used", n_used);

    otherwise   % "nn"
        if isempty(opts.X_s) || isempty(opts.X_t)
            error("TEN:I_ALIGNEMBED:MissingExpression", ...
                "solver=""nn"" needs the expression matrices: pass X_s and " + ...
                "X_t. The network maps each gene's profile across cells to " + ...
                "its coordinates, so unlike the spectral solver it cannot " + ...
                "work from W alone.");
        end
        [P_s, P_t, info] = ten.i_nnembed(W, ng, n_dim, opts.X_s, opts.X_t, ...
            n_steps=opts.n_steps, lr=opts.lr, seed=opts.seed, ...
            gradient=opts.gradient, degree=opts.degree, verbose=opts.verbose);
        info.solver = "nn";
end

end
