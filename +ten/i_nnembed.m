function [P_s, P_t, info] = i_nnembed(W, ng, n_dim, X_s, X_t, opts)
%I_NNEMBED  Neural-network manifold alignment of a two-block graph.
%   [P_S, P_T] = TEN.I_NNEMBED(W, NG, N_DIM, X_S, X_T) embeds the genes of a
%   2*NG-by-2*NG block alignment graph by training two feedforward networks to
%   minimise the graph-Laplacian loss, mirroring the published scTenifoldXct
%   training loop.
%
%   This is NOT an alternative solver for the same problem TEN.I_LAPLACIANEMBED
%   solves. Both minimise trace(P'*L*P), but over different parameter spaces:
%
%     TEN.I_LAPLACIANEMBED - P is free; the constrained minimum is attained in
%                            closed form by the bottom eigenvectors of L. This
%                            is the global optimum.
%     TEN.I_NNEMBED        - P is the Stiefel retraction of the outputs of two
%                            MLPs applied to the expression matrices, and the
%                            network weights are optimised. P is confined to the
%                            image of that function class, so the attained loss
%                            is generally higher.
%
%   The restriction acts as a prior: because a shared function generates the
%   coordinates, genes with similar expression profiles receive similar
%   embeddings. The spectral route sees expression only through the thresholded
%   GRN encoded in W. Disagreement between the two is therefore informative
%   about that prior rather than evidence of a defect in either.
%
%   INPUTS:
%     W     - 2*NG-by-2*NG symmetric block weight matrix
%     ng    - number of genes (half the side length of W)
%     n_dim - embedding dimensions
%     X_s   - NG-by-nCellsSource expression matrix (rows are genes)
%     X_t   - NG-by-nCellsTarget expression matrix (rows are genes)
%
%   NAME-VALUE ARGUMENTS:
%     n_steps    - Adam iterations (default 1000)
%     lr         - Adam learning rate (default 0.01)
%     seed       - RNG seed for weight initialisation (default 0). Exposed so
%                  that seed sensitivity can be measured; the published loop
%                  fixes it at 0.
%     gradient   - "riemannian" (default) or "euclidean". See below; prefer
%                  the default, and use "euclidean" only to reproduce runs
%                  recorded before this option existed.
%     degree     - "signed" (default) or "abs". Which degree the Laplacian
%                  uses; see DEGREE CONVENTION below.
%     nsIter     - Newton-Schulz iterations for the differentiable Stiefel
%                  retraction, "euclidean" path only (default 10)
%     verbose    - print training progress (default false)
%
%   OUTPUTS:
%     P_s  - NG-by-N_DIM source gene coordinates
%     P_t  - NG-by-N_DIM target gene coordinates
%     info - struct with diagnostics:
%              .loss          final loss value
%              .lossTrace     loss at every step
%              .orthErr       max(abs(P'*P - I)) for the returned P
%              .ncomp         connected components of the alignment graph
%              .nullfrac      fraction of the embedding's energy lying in the
%                             Laplacian's null space; ~1 means no structure
%              .nulldominant  true when ncomp > n_dim, i.e. a zero-loss
%                             null-space embedding is guaranteed to exist
%              .degenerate    true when the embedding is measurably null space
%
%   GRADIENT PATHS
%
%   "riemannian" (default) follows the published implementation. Each step
%   retracts the network outputs with an exact economy SVD, forms the analytic
%   Euclidean gradient of the loss, 2*L*P/3000, projects it onto the Stiefel
%   tangent space, and backpropagates that Riemannian gradient into the network
%   weights.
%
%   The reference does this with proj_outputs.backward(rgrad), which carries
%   the cotangent through torch.svd. MATLAB's svd does not accept a dlarray, so
%   the retraction is computed outside the autodiff graph and its chain rule
%   supplied analytically: differentiating P = U*V' gives
%       dP = U*K*V' + (I - U*U')*dY*V*inv(S)*V',
%       K_ij = (A_ij - A_ji)/(s_i + s_j),   A = U'*dY*V,
%   in which the usual 1/(s_j^2 - s_i^2) poles cancel, so repeated singular
%   values are harmless. The adjoint of that map is applied to rgrad before the
%   network backward pass, making this equivalent to the reference rather than
%   an approximation of it. Earlier versions seeded the backward pass at the
%   pre-retraction outputs instead, which was the residual MATLAB-Python gap.
%
%   "euclidean" is the original path, kept only for reproducing earlier runs.
%   It differentiates the loss straight through a fixed-step Newton-Schulz
%   retraction. That is measurably worse, and not by a small margin. On a
%   120-gene synthetic fixture the published loop reaches loss 0.0895 while
%   this path reaches 0.15 to 0.20. The reason is that the network outputs are
%   severely ill-conditioned - a representative normalized spectrum runs from
%   1.0 down to 2.7e-4 - so ten Newton-Schulz steps leave P'*P off the identity
%   by ~0.97, and the polar factor's derivative is badly behaved besides.
%   Raising nsIter does not rescue it: 30 iterations orthogonalize properly but
%   give a worse loss still, because the defect is the differentiated
%   retraction rather than the retraction's accuracy.
%
%   info.orthErr is measured on the returned embedding, which is always
%   retracted by an exact SVD, so it is near machine precision on both paths
%   and says nothing about training-time orthogonality. Do not read it as
%   evidence that the "euclidean" path is well behaved.
%
%   DEGREE CONVENTION
%
%   "signed" (default) builds L = diag(sum(W)) - W, reproducing the reference's
%   scipy.sparse.csgraph.laplacian(w, normed=False). "abs" uses sum(abs(W)),
%   which was the previous MATLAB behaviour.
%
%   The two agree exactly on a non-negative W and differ otherwise. This
%   matters because a pcrnet GRN is signed: under "signed" the Laplacian is
%   then indefinite, the loss can go negative, and the null space is no longer
%   spanned by the connected-component indicators, which is what the degeneracy
%   diagnostics below assume. The reference does not hit this because core.py
%   adds 1 to every element of each GRN block before assembly, making W
%   non-negative and the two conventions identical. MATLAB does not yet apply
%   that offset, so on a signed W the "signed" convention is faithful to the
%   reference's formula but not to the regime the reference actually runs in.
%   Use "abs" if you need a PSD Laplacian and guaranteed-meaningful degeneracy
%   diagnostics on a signed W.
%
%   RETURNED EMBEDDING
%
%   The reference returns the projection computed inside the final iteration,
%   before that iteration's optimiser step, so its embedding reflects
%   n_steps-1 updates. That is reproduced here on the "riemannian" path.
%
%   DEGENERACY
%
%   On a disconnected graph the loss is minimised on null-space vectors and
%   training can converge there. Unlike the spectral route, which sizes its
%   eigenvalue request past the null space and so cannot return one, this
%   function will happily return a null-space embedding: nonzero, orthonormal,
%   entirely ordinary-looking coordinates carrying no information about W.
%
%   The condition is therefore measured rather than inferred from the loss. The
%   null space of L is spanned by the indicator vectors of the connected
%   components, so info.nullfrac is the fraction of the embedding's energy
%   lying in that span, computed exactly. A loss threshold cannot do this job:
%   the loss has no intrinsic scale, and an embedding can sit almost wholly in
%   the null space at a loss orders of magnitude above any fixed tolerance. On
%   the glyco validation fixture, 168 components at n_dim=20 gave loss 4.1e-4 -
%   far above the 1e-8 tolerance this function used to test against, and far
%   too small to mean anything.
%
%   info.nulldominant is the structural companion: when ncomp > n_dim an
%   exactly-zero-loss embedding exists inside the null space, so the objective
%   cannot distinguish structure from noise no matter how well it is optimised.
%   That is a property of the graph alone and needs no tolerance.
%
%   Reading info.nullfrac: it is bounded above by min(ncomp, n_dim)/n_dim, so a
%   small value is not automatically reassuring and a moderate one is not
%   automatically alarming. On a connected graph the null space is the single
%   constant vector, which the embedding will normally contain, giving
%   nullfrac ~ 1/n_dim as the healthy baseline - 0.2 at n_dim=5, 0.05 at
%   n_dim=20. Only values approaching 1 indicate an embedding with nothing in
%   it, and reaching 1 requires ncomp >= n_dim in the first place.
%
%   Requires Deep Learning Toolbox (R2022b+).
%
% see also: TEN.I_ALIGNEMBED, TEN.I_LAPLACIANEMBED, TEN.XCT.XCTMAIN_NN

arguments
    W {mustBeNumeric}
    ng (1, 1) double
    n_dim (1, 1) double
    X_s {mustBeNumeric}
    X_t {mustBeNumeric}
    opts.n_steps (1, 1) double = 1000
    opts.lr (1, 1) double = 0.01
    opts.seed (1, 1) double = 0
    opts.gradient (1, 1) string ...
        {mustBeMember(opts.gradient, ["riemannian", "euclidean"])} = "riemannian"
    opts.degree (1, 1) string ...
        {mustBeMember(opts.degree, ["signed", "abs"])} = "signed"
    opts.nsIter (1, 1) double = 10
    opts.verbose (1, 1) logical = false
end

if ~(exist('dlarray', 'builtin') || exist('dlarray', 'file'))
    error("TEN:I_NNEMBED:NoDeepLearning", ...
        "Neural-network alignment requires the Deep Learning Toolbox " + ...
        "(R2022b+). Use solver=""spectral"" for a toolbox-free alternative.");
end
if size(X_s, 1) ~= ng || size(X_t, 1) ~= ng
    error("TEN:I_NNEMBED:BadSizes", ...
        "X_S and X_T must both have %d rows (one per gene).", ng);
end

n_steps = round(opts.n_steps);

% ── Unnormalised Laplacian, dense single for the autodiff path ───────────
% "signed" reproduces scipy.sparse.csgraph.laplacian(w, normed=False), which is
% what the reference calls. "abs" is the previous MATLAB convention; see the
% DEGREE CONVENTION section of the help.
if opts.degree == "abs"
    d = full(sum(abs(W), 2));
else
    d = full(sum(W, 2));
end
L_dense = single(full(spdiags(d, 0, 2*ng, 2*ng) - W));

% ── Initialise both networks ─────────────────────────────────────────────
rngState = rng;                 % restore the caller's stream on exit
cleanup = onCleanup(@() rng(rngState));
rng(opts.seed);

ps = i_initparams(size(X_s, 2), n_dim);
pt = i_initparams(size(X_t, 2), n_dim);

% dlarray rejects sparse input with an error that does not name the culprit,
% so densify here. The callers already pass full matrices; this is for direct
% callers.
Xs_dl = dlarray(single(full(X_s)));
Xt_dl = dlarray(single(full(X_t)));
L_dl = dlarray(L_dense);

% ── Training loop ────────────────────────────────────────────────────────
[avg_s, avgSq_s] = deal([]);
[avg_t, avgSq_t] = deal([]);
lossTrace = zeros(n_steps, 1);

useRiemannian = opts.gradient == "riemannian";
P = [];

for step = 1:n_steps
    if useRiemannian
        [lossVal, grads, P] = dlfeval(@i_gradfnriem, ps, pt, Xs_dl, Xt_dl, L_dense);
    else
        [lossVal, grads] = dlfeval(@i_gradfn, ps, pt, Xs_dl, Xt_dl, L_dl, opts.nsIter);
    end
    [ps, avg_s, avgSq_s] = adamupdate(ps, grads{1}, avg_s, avgSq_s, step, opts.lr);
    [pt, avg_t, avgSq_t] = adamupdate(pt, grads{2}, avg_t, avgSq_t, step, opts.lr);

    lossTrace(step) = double(extractdata(lossVal));
    if opts.verbose && (step == 1 || mod(step, 100) == 0)
        fprintf('[i_nnembed]   step %4d / %d   loss = %.6f\n', ...
            step, n_steps, lossTrace(step));
    end
end

% ── Final embedding ──────────────────────────────────────────────────────
% The reference returns the projection computed inside the last iteration,
% before that iteration's optimiser step, so the returned embedding reflects
% n_steps-1 updates rather than n_steps. That is reproduced here. The
% "euclidean" path has no exactly-retracted in-loop P to return - its
% Newton-Schulz iterate is far off the manifold - so it keeps the post-hoc
% exact retraction it always used.
if ~useRiemannian || isempty(P)
    outputs = [extractdata(i_netfwd(Xs_dl, ps)); extractdata(i_netfwd(Xt_dl, pt))];
    [U, ~, V] = svd(outputs, "econ");
    P = U*V';
end
P = double(P);

P_s = P(1:ng, :);
P_t = P(ng+1:end, :);

% ── Diagnostics ──────────────────────────────────────────────────────────
gram = P'*P;
info.orthErr = max(abs(gram - eye(size(gram))), [], "all");
info.loss = sum(P.*(double(L_dense)*P), "all")/3000;
info.lossTrace = lossTrace;

% Measure how much of the embedding lies in the null space, rather than
% inferring it from the loss. See the DEGENERACY section of the help above for
% why a loss tolerance cannot detect this.
[info.ncomp, comp] = i_ncomponents(W);
info.nullfrac = i_nullfraction(P, comp, info.ncomp);
info.nulldominant = info.ncomp > n_dim;
info.degenerate = info.ncomp > 1 && info.nullfrac > 0.99;

advice = "Lower the GRN edge threshold or filter low-expression genes to " + ...
    "connect the graph, or raise n_dim above the component count.";
if info.degenerate
    warning("TEN:I_NNEMBED:Degenerate", ...
        "The alignment graph has %d connected components and %.1f%% of the " + ...
        "embedding's energy lies in the Laplacian's null space, so these " + ...
        "coordinates carry essentially no information about W. " + advice, ...
        info.ncomp, 100*info.nullfrac);
elseif info.nulldominant
    warning("TEN:I_NNEMBED:NullDominant", ...
        "The alignment graph has %d connected components but n_dim is only " + ...
        "%d, so an exactly-zero-loss embedding exists entirely within the " + ...
        "null space and the objective cannot distinguish it from structure. " + ...
        "Measured null-space energy: %.1f%%. " + advice, ...
        info.ncomp, n_dim, 100*info.nullfrac);
end

end


%% ---- Riemannian path: exact retraction, projected gradient ----
function [loss, grads, P] = i_gradfnriem(ps, pt, Xs_dl, Xt_dl, L)
% Mirrors ManifoldAlignmentNet.train: retract exactly, take the Euclidean
% gradient at the retracted point, project it onto the Stiefel tangent space,
% and backpropagate that Riemannian gradient through the retraction into the
% network weights.
%
% The retraction is computed outside the autodiff graph because MATLAB's svd
% does not accept a dlarray. The chain rule through it is instead supplied
% analytically by i_polarback, so this is equivalent to the reference's
% proj_outputs.backward(rgrad), which does differentiate through torch.svd.

outputs = [i_netfwd(Xs_dl, ps); i_netfwd(Xt_dl, pt)];

Y = double(extractdata(outputs));
[U, S, V] = svd(Y, "econ");
P = U*V';

% Euclidean gradient of trace(P'*L*P)/3000, analytic since L is symmetric.
LP = double(L)*P;
egrad = (2/3000).*LP;

% Tangent-space projection, as in stiefel.proj_stiefel. The second term is
% written as Z - P*(P'*Z) rather than (I - P*P')*Z so that no 2*ng-by-2*ng
% identity is ever formed.
skewPart = 0.5*(P'*egrad - egrad'*P);
rgrad = P*skewPart + (egrad - P*(P'*egrad));

% Pull the cotangent back through the retraction, then seed the network
% backward pass with it.
Ybar = i_polarback(U, S, V, rgrad);

surrogate = sum(outputs.*dlarray(single(Ybar)), "all");
grads = dlgradient(surrogate, {ps, pt});

loss = dlarray(single(sum(P.*LP, "all")/3000));

end


%% ---- adjoint of the polar factor P = U*V' of Y = U*S*V' ----
function Ybar = i_polarback(U, S, V, G)
% Given the cotangent G = df/dP, return df/dY.
%
% Differentiating P = U*V' gives
%     dP = U*K*V' + (I - U*U')*dY*V*inv(S)*V',
%     K_ij = (A_ij - A_ji)/(s_i + s_j),   A = U'*dY*V,
% where the usual 1/(s_j^2 - s_i^2) poles cancel against the (s_j - s_i)
% factor, so no eigenvalue-gap term survives and repeated singular values are
% harmless. Taking the adjoint of that map gives the expression below.
%
% Verified against central differences to ~1e-9 relative, for both generic and
% Stiefel-tangent cotangents.

s = diag(S);
tol = max(s)*eps(class(s))*numel(s);
s = max(s, tol);                % guard a rank-deficient Y

B = U'*G*V;
C = B./(s + s');                % elementwise, outer sum of singular values

Ybar = U*(C - C')*V' + (G - U*(U'*G))*V*diag(1./s)*V';

end


%% ---- Euclidean path: differentiate through the retraction ----
function [loss, grads] = i_gradfn(ps, pt, Xs_dl, Xt_dl, L_dl, nsIter)
outputs = [i_netfwd(Xs_dl, ps); i_netfwd(Xt_dl, pt)];

% Differentiable Stiefel retraction. dlsvd does not exist, so the polar factor
% is approached by Newton-Schulz, which uses only matmul and scalar ops.
P = outputs ./ sqrt(sum(outputs.^2, "all"));
for k = 1:nsIter
    P = 1.5.*P - 0.5.*(P*(P'*P));
end

loss = sum(P.*(L_dl*P), "all")/3000;
grads = dlgradient(loss, {ps, pt});
end


%% ---- forward pass through one 3-layer network ----
function out = i_netfwd(X, p)
h1 = sigmoid(X*p.W1 + p.b1);
h2 = sigmoid(h1*p.W2 + p.b2);
out = h2*p.W3 + p.b3;
end


%% ---- parameters: n_cells -> 4*n_h -> n_h -> n_dim ----
function p = i_initparams(n_cells, n_dim)
% Layer widths follow create_models: n_h = gmean([n_cells, n_dim]) truncated to
% an integer, which for two values is sqrt(n_cells*n_dim), and a = 4.

n_h = max(floor(sqrt(n_cells*n_dim)), 1);
H1 = 4*n_h;
H2 = n_h;

[p.W1, p.b1] = i_linearinit(n_cells, H1);
[p.W2, p.b2] = i_linearinit(H1, H2);
[p.W3, p.b3] = i_linearinit(H2, n_dim);
end


function [W, b] = i_linearinit(fan_in, fan_out)
% torch.nn.Linear's default reset_parameters: kaiming_uniform_(a=sqrt(5)) on
% the weight, which reduces to U(-1/sqrt(fan_in), 1/sqrt(fan_in)), and the bias
% drawn from the same bound. Note the reference relies on this default rather
% than initialising explicitly.

bound = 1/sqrt(fan_in);
W = dlarray((rand(fan_in, fan_out, "single")*2 - 1).*bound);
b = dlarray((rand(1, fan_out, "single")*2 - 1).*bound);
end


%% ---- connected components of the graph underlying W ----
function [ncomp, comp] = i_ncomponents(W)
A = abs(W) > 0;
A = A | A';
try
    comp = conncomp(graph(A, "omitselfloops"));
    ncomp = max(comp);
catch
    comp = ones(1, size(W, 1));
    ncomp = 1;
end
end


%% ---- fraction of the embedding's energy in the Laplacian's null space ----
function frac = i_nullfraction(P, comp, ncomp)
% The null space of an unnormalised Laplacian is spanned by the indicator
% vectors of the connected components. Projecting onto their normalised span
% is a per-component sum, so no basis is formed explicitly.

total = sum(P.^2, "all");
if total <= 0
    frac = 0;
    return;
end

n = size(P, 1);
M = sparse(comp(:), (1:n)', 1, ncomp, n);
compSum = M*P;                  % ncomp-by-n_dim component sums
compSize = full(sum(M, 2));     % ncomp-by-1
frac = sum((compSum.^2)./compSize, "all")/total;

end
