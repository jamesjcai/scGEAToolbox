function [A, Tedges, mps, ordidx] = tngrn(X, genelist, opts)
% TNGRN  Tensor-network (MPS Born machine) GRN inference for a gene panel.
%
%   A = ten.tngrn(X)
%   A = ten.tngrn(X, genelist)
%   [A, Tedges, mps, ordidx] = ten.tngrn(X, genelist, Name, Value, ...)
%
%   Panel-scale companion to sctenifoldnet / sc_grn backends. Where pcrnet gives
%   a pairwise regression network genome-wide, ten.tngrn takes a curated module
%   (a TF set / pathway, typically <= ~20 genes), learns a matrix-product-state
%   Born machine of the joint activity distribution (i_mpsborntrain), then reads
%   edges off it as mutual information of the model marginals (qtm.mpsqmi;
%   classical MI, not von Neumann QMI -- see that file) with a permutation
%   null. Beyond pairwise edges it exposes (i) virtual knockout by site clamping
%   -- scTenifoldKnk logic as an O(1) tensor op -- and (ii, TODO) triadic
%   interaction information, neither of which the regression network can give.
%
%   X          genes x cells expression (log-normalised; PFlog1pPF fine), panel-sized
%   genelist   string/cellstr of gene names (default "g1".."gN")
%
%   Name-Value
%     LocalDim   d, discretization levels per gene            (default 2)
%     BondDim    Dmax, MPS bond-dimension cap                 (default 20)
%     Order      "hilbert" | "none" | "input" | index vector  (default "hilbert")
%     Binarize   "gmm" | "quantile"                           (default "gmm")
%     NumSweeps  DMRG sweeps                                  (default 12)
%     LearnRate  gradient step size                           (default 0.1)
%     NumGrad    gradient steps per two-site tensor           (default 3)
%     NumPerm    permutations for the QMI null                (default 1000)
%     Alpha      significance level (two one-tailed tests)    (default 0.10)
%     Knockout   gene name/index to clamp to inactive, or []  (default [])
%     EdgesOnly  return only the QMI adjacency A (for sc_grn) (default false)
%
%   Outputs
%     A        [N x N] signed MI adjacency: +MI for right-tail (association),
%              -MI for left-tail (mutual exclusivity / repression), 0 if n.s.
%     Tedges   table Source,Target,QMI_bits,pRight,pLeft,Sign
%              (QMI_bits column holds classical MI bits; name kept for API stability)
%     mps      trained MPS (cell of [Dl x d x Dr] cores), gene order = ordidx
%     ordidx   chain ordering applied to genes
%
%   Estimator note: the point MI is read from the trained MPS marginals
%   (denoised, higher-order-aware); the null is empirical (shuffled d-level
%   marginals, same estimator family as qtm.BinPairMI). On ~10^4 cells the
%   two-site marginals are well sampled, so it is a calibration null, not the
%   point estimator. (MI here is classical Shannon MI, not von Neumann QMI --
%   see qtm.mpsqmi.)
%
%   See also SCTENIFOLDNET, SCTENIFOLDKNK, I_MPSBORNTRAIN, qtm.MPSQMI, SC_GRN.

arguments
    X {mustBeNumeric}
    genelist string = string.empty
    opts.LocalDim  (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.LocalDim,2)} = 2
    opts.BondDim   (1,1) double {mustBeInteger} = 20
    opts.Order     = "hilbert"
    opts.Binarize  (1,1) string = "gmm"
    opts.NumSweeps (1,1) double = 12
    opts.LearnRate (1,1) double = 0.1
    opts.NumGrad   (1,1) double = 3
    opts.NumPerm   (1,1) double = 1000
    opts.Alpha     (1,1) double = 0.10
    opts.Knockout  = []
    opts.EdgesOnly (1,1) logical = false
end

[G, ~] = size(X);
if isempty(genelist), genelist = "g" + string(1:G); end
genelist = string(genelist(:));
assert(numel(genelist) == G, 'genelist length must match rows of X.');
if G > 25
    warning('ten:tngrn:panelSize', ...
        ['N=%d genes: tngrn uses an exact d^N joint readout and is panel-scale. ' ...
         'Restrict to a module, or switch qtm.mpsqmi to the RDM route.'], G);
end

d = opts.LocalDim;

% ---- 1. discretize genes to codes 1..d ---------------------------------
V = i_discretize(X, d, opts.Binarize);          % [cells x genes]

% ---- 2. chain ordering (MPS is 1D; correlated genes adjacent) ----------
ordidx = i_chain_order(V, opts.Order);
V = V(:, ordidx);
genes = genelist(ordidx);

% ---- 3. train MPS Born machine (two-site DMRG) -------------------------
mps = ten.i_mpsborntrain(V, d, opts.BondDim, ...
        struct('nsweeps', opts.NumSweeps, 'lr', opts.LearnRate, ...
               'ngrad', opts.NumGrad, 'tol', 1e-8));

% ---- optional virtual knockout: clamp a site to inactive ---------------
if ~isempty(opts.Knockout)
    ko  = i_resolve_gene(opts.Knockout, genes);
    mps = i_mps_clamp(mps, ko, 1);              % qtm.mpsqmi renormalises => conditional
end

% ---- 4. point QMI from the trained MPS ---------------------------------
Q = qtm.mpsqmi(mps, d);                          % [N x N] bits

if opts.EdgesOnly
    A = Q; Tedges = table(); return
end

% ---- 5. empirical permutation null -------------------------------------
[pR, pL] = i_permnull(V, Q, d, opts.NumPerm);

% ---- 6. signed thresholded adjacency + edge table ----------------------
N   = numel(genes);
sig = (pR < opts.Alpha) | (pL < opts.Alpha);
sgn = double(pR <= pL) - double(pL < pR);        % + assoc / - exclusivity
A   = zeros(N);
A(sig) = Q(sig) .* sgn(sig);
A = (A + A.') / 2;                               % symmetric (direction needs clamping)

[si, ti] = find(triu(sig, 1));
lin = sub2ind([N N], si, ti);
Tedges = table(genes(si), genes(ti), Q(lin), pR(lin), pL(lin), sgn(lin), ...
    'VariableNames', {'Source', 'Target', 'QMI_bits', 'pRight', 'pLeft', 'Sign'});
end % tngrn


% ========================================================================
% Local helpers
% ========================================================================

function V = i_discretize(X, d, method)
% Per-gene discretization to codes 1..d. GMM branch mirrors Sanz Larrarte:
% fit a d-component GMM per gene, order components by mean (code 1 = inactive),
% assign by max posterior.
[G, m] = size(X);
V = ones(m, G);
X = double(full(X));
for g = 1:G
    x = X(g, :).';
    switch lower(method)
        case "gmm"
            gm = fitgmdist(x, d, 'Replicates', 3, 'RegularizationValue', 1e-6);
            [~, ord] = sort(gm.mu);              % ord(r) = component with r-th smallest mean
            [~, lab] = max(posterior(gm, x), [], 2);
            remap(ord) = 1:d;                    %#ok<AGROW>
            V(:, g) = remap(lab);
        case "quantile"
            edges = quantile(x, linspace(0, 1, d + 1));
            V(:, g) = min(max(discretize(x, edges), 1), d);
    end
end
end

function ordidx = i_chain_order(V, order)
% Hilbert order: PCA genes to 2D, map to 1D via a Hilbert curve so correlated
% genes stay adjacent (keeps bond entropy low). Alternatives: minimise summed
% bond entropy, or move to a tree tensor network.
N = size(V, 2);
if isnumeric(order), ordidx = order(:).'; return; end
switch lower(string(order))
    case "hilbert"
        [~, sc] = pca(double(V).', 'NumComponents', 2);   % genes x 2
        ordidx = i_hilbert_sort(sc);
    otherwise
        ordidx = 1:N;                            % "none" / "input"
end
end

function idx = i_hilbert_sort(xy)
% PLACEHOLDER: rank 2D points along a Hilbert curve. Replace with a true
% d2xy/xy2d bit-interleave on a quantised grid (a dozen lines). sortrows keeps
% the pipeline runnable meanwhile.
[~, idx] = sortrows(xy);
idx = idx(:).';
end

function ko = i_resolve_gene(g, genes)
if isnumeric(g), ko = g; else, ko = find(genes == string(g), 1); end
assert(~isempty(ko), 'Knockout gene not found in panel.');
end

function mps = i_mps_clamp(mps, site, code)
% Virtual knockout: project a site onto a fixed physical code. Renormalisation
% is deferred to qtm.mpsqmi (which divides by sum(P)), yielding the conditional.
A = mps{site};
keepvec = zeros(1, size(A, 2)); keepvec(code) = 1;
mps{site} = A .* reshape(keepvec, [1 numel(keepvec) 1]);
end

function [pR, pL] = i_permnull(V, Q, d, nperm)
% Empirical null: permute each gene column independently to destroy dependence,
% recompute pairwise MI over fixed d-level bins (same estimator as qtm.mpsqmi's
% marginals, so model-vs-empirical joint is the only difference), and count
% exceedances. Parfor-able over b for large panels.
[m, N] = size(V);
gR = zeros(N); gL = zeros(N);
for b = 1:nperm
    Vp = V;
    for g = 1:N, Vp(:, g) = V(randperm(m), g); end
    for i = 1:N
        for j = i + 1:N
            qn = i_mi2(Vp(:, i), Vp(:, j), d);
            gR(i, j) = gR(i, j) + (qn >= Q(i, j));
            gL(i, j) = gL(i, j) + (qn <= Q(i, j));
        end
    end
end
pR = (gR + gR.' + 1) / (nperm + 1);
pL = (gL + gL.' + 1) / (nperm + 1);
pR(1:N + 1:end) = 1; pL(1:N + 1:end) = 1;
end

function mi = i_mi2(a, b, d)
% Shannon MI (bits) of two integer-code vectors over fixed d levels.
Pij = accumarray([a b], 1, [d d]); Pij = Pij / sum(Pij(:));
Pi = sum(Pij, 2); Pj = sum(Pij, 1);
ent = @(pv) -sum(pv(pv > 0) .* log2(pv(pv > 0)));
mi = ent(Pi) + ent(Pj(:)) - ent(Pij(:));
end
