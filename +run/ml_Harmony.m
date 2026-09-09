function [Z_corr, R, model] = ml_Harmony(Z, batch, opts)
%ML_HARMONY Integrate batches by correcting a low-dimensional embedding.
%   Z_corr = ML_HARMONY(Z, batch) removes batch effects from an embedding
%   (usually PCA scores) while preserving the biological structure, by
%   alternating soft clustering with a per-cluster, per-batch linear
%   correction. This is a native MATLAB implementation of Harmony
%   (Korsunsky et al., Nat Methods 2019, PMID:31740819).
%
%   The premise is that a batch effect is a shift that a cell type shares
%   with every other cell of its batch, so it can be estimated cluster by
%   cluster and subtracted. The clustering is deliberately soft and is
%   penalised for being batch-pure: a cluster made of one batch alone gives
%   no evidence about what the shift is, and the diversity penalty pushes
%   the assignment away from that degenerate solution.
%
%   INPUTS:
%     Z     - numCells-by-numDims embedding. Rows are cells.
%     batch - numCells-by-1 batch labels (numeric, string, cellstr or
%             categorical).
%
%   NAME-VALUE:
%     K        - number of soft clusters (default min(100, round(n/30)),
%                at least 2).
%     Sigma    - soft assignment bandwidth on the cosine distance. Default
%                [] picks it from the data so that a cell spreads over
%                PERPLEXITY clusters on average; pass a number to fix it.
%                The published default of 0.1 assumes cosine distances on
%                scaled PCs and produces near-hard assignments on anything
%                else, which stalls the correction (see note below).
%     Perplexity - target effective number of clusters per cell used to
%                pick Sigma automatically (default min(K, 8)). Ignored when
%                Sigma is given.
%     Theta    - diversity penalty strength (default 2). 0 disables the
%                penalty and reduces the method to per-cluster centering;
%                larger values force harder mixing across batches.
%     Lambda   - ridge penalty on the batch coefficients (default 1). The
%                intercept is never penalised and never removed.
%     MaxIter  - maximum outer (cluster + correct) iterations (default 10).
%     Tol      - relative change of Z_corr at which to stop (default 1e-4).
%     Verbose  - print progress (default false).
%     Seed     - RNG seed for the k-means initialisation (default [], no
%                reseed).
%
%   OUTPUTS:
%     Z_corr - numCells-by-numDims corrected embedding, in the same units
%              as Z.
%     R      - numCells-by-K soft cluster responsibilities, rows sum to 1.
%     model  - struct with the cluster centroids (Y), the observed and
%              expected batch counts per cluster (O, E), the batch levels
%              and the parameters used.
%
%   EXAMPLE:
%     [~, sc] = pca(log1p(sc_norm(sce.X))', "NumComponents", 50);
%     sc = run.ml_Harmony(sc, sce.c_batch_id);
%     sce.s = sc_umap(sc');
%
%   Note on LAMBDA: before 26.4 this argument named the diversity penalty
%   exponent. It is now the ridge penalty, matching the reference
%   implementation; use THETA for the diversity penalty.
%
%   Note on SIGMA: the correction is only identifiable while clusters
%   contain more than one batch. If SIGMA is small relative to the spread
%   of the cosine distances, the responsibilities collapse to a hard
%   assignment, every cluster becomes batch-pure, and the batch term of the
%   regression goes to zero -- the function then returns Z unchanged. The
%   automatic bandwidth exists to keep the assignment soft enough that this
%   cannot happen silently.
%
%   See also SC_UMAP, SC_TSNE, GUI.CALLBACK_HARMONY.

arguments
    Z (:,:) {mustBeNumeric, mustBeNonempty}
    batch
    opts.K (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    opts.Sigma double {mustBeScalarOrEmpty, mustBePositive} = []
    opts.Perplexity (1,1) double {mustBeNonnegative} = 0
    opts.Theta (1,1) double {mustBeNonnegative} = 2
    opts.Lambda (1,1) double {mustBeNonnegative} = 1
    opts.MaxIter (1,1) double {mustBePositive, mustBeInteger} = 10
    opts.Tol (1,1) double {mustBePositive} = 1e-4
    opts.Verbose (1,1) logical = false
    opts.Seed = []
end

[n, d] = size(Z);
if numel(batch) ~= n
    error("ml_Harmony:batchSize", ...
        "BATCH has %d entries but Z has %d rows. Supply one batch label per cell.", ...
        numel(batch), n);
end

if ~iscategorical(batch)
    batch = categorical(batch(:));
end
[batchLevels, ~, lbl] = unique(batch(:));
numBatches = numel(batchLevels);
if numBatches < 2
    warning("ml_Harmony:singleBatch", ...
        "All cells share one batch label; Z is returned unchanged.");
    Z_corr = double(Z);
    R = ones(n, 1);
    model = struct("Y", [], "O", [], "E", [], ...
        "batchLevels", batchLevels, "params", opts);
    return;
end

K = opts.K;
if K == 0
    K = min(100, max(2, round(n/30)));
end
K = min(K, n);

if ~isempty(opts.Seed)
    rng(opts.Seed);
end

% Work in dims-by-cells throughout: the correction is a small regression per
% cluster and this orientation keeps those solves contiguous.
Zt = double(Z).';
Zcorr = Zt;
Zcos = i_normalizecols(Zcorr);

% One-hot batch design, plus an intercept row for the correction step.
phi = false(numBatches, n);
phi(sub2ind([numBatches, n], lbl(:).', 1:n)) = true;
phi = double(phi);
phiMoe = [ones(1, n); phi];
ridge = diag([0, opts.Lambda*ones(1, numBatches)]);
priorBatch = sum(phi, 2) / n;

% k-means on the cosine-normalised data gives centroids on the unit sphere,
% which is the space the responsibilities are computed in.
[~, Y] = kmeans(Zcos.', K, "MaxIter", 100, "Replicates", 3, ...
    "Display", "off");
Y = i_normalizecols(Y.');

% Pick the bandwidth from the data unless the caller fixed it. A cell that
% belongs to one cluster alone carries no information about how its batch
% differs from the others, so the assignment has to stay soft.
sigma = opts.Sigma;
if isempty(sigma)
    perplexity = opts.Perplexity;
    if perplexity == 0
        perplexity = min(K, 8);
    end
    sigma = i_autosigma(2*(1 - Y.'*Zcos), min(perplexity, K));
    if opts.Verbose
        fprintf("ml_Harmony: Sigma = %.4g (target perplexity %.3g)\n", ...
            sigma, min(perplexity, K));
    end
end

R = i_softassign(2*(1 - Y.'*Zcos), sigma);
O = R*phi.';                    % K-by-numBatches observed batch mass
E = sum(R, 2)*priorBatch.';     % expected mass under batch independence

for iter = 1:opts.MaxIter
    [R, Y, O, E] = i_cluster(Zcos, R, O, E, phi, priorBatch, sigma, ...
        opts.Theta);

    % ---- Correction: per-cluster weighted ridge, intercept retained ----
    % Every cluster proposes a shift for each batch; a cell is moved by the
    % sum of those shifts weighted by how much it belongs to each cluster.
    % Starting from Zt rather than the running estimate keeps the total
    % correction a single linear map of the original embedding.
    Znew = Zt;
    live = find(sum(R, 2) > 1e-8);
    for k = live(:).'
        phiRk = phiMoe.*R(k, :);
        W = (phiRk*phiMoe.' + ridge) \ (phiRk*Zt.');
        W(1, :) = 0;                    % keep the intercept: it is biology
        Znew = Znew - W.'*phiRk;
    end

    rel = norm(Znew - Zcorr, "fro")/(norm(Zcorr, "fro") + eps);
    Zcorr = Znew;
    Zcos = i_normalizecols(Zcorr);
    if opts.Verbose
        fprintf("ml_Harmony iter %2d: relative change %.3e\n", iter, rel);
    end
    if rel < opts.Tol
        break;
    end
end

Z_corr = Zcorr.';
R = R.';
model = struct("Y", Y, "O", O, "E", E, ...
    "batchLevels", batchLevels, "params", opts, "Sigma", sigma);
if size(Z_corr, 2) ~= d
    error("ml_Harmony:internal", "Corrected embedding lost its dimensions.");
end

end


function [R, Y, O, E] = i_cluster(Zcos, R, O, E, phi, priorBatch, sigma, theta)
% Soft k-means with the diversity penalty, updated over random blocks of
% cells. Two details keep this stable, and the method does nothing without
% either of them. The penalty is smoothed as ((E+1)/(O+1))^theta, so that a
% cluster momentarily empty of one batch does not exert an unbounded pull.
% And each block's mass is taken out of O and E before its responsibilities
% are recomputed, then added back; updating every cell at once against one
% frozen O makes the penalty self-reinforcing, and the clusters collapse
% onto a handful of centroids that are each batch-pure -- at which point
% the correction is exactly zero and the whole method silently no-ops.

n = size(Zcos, 2);
numBlocks = 5;
maxIter = 20;

for iter = 1:maxIter
    Y = i_normalizecols(Zcos*R.');
    dist = 2*(1 - Y.'*Zcos);
    base = exp(-(dist - min(dist, [], 1))/sigma);

    Rprev = R;
    order = randperm(n);
    edges = round(linspace(0, n, numBlocks+1));
    for b = 1:numBlocks
        idx = order(edges(b)+1:edges(b+1));
        O = O - R(:, idx)*phi(:, idx).';
        E = E - sum(R(:, idx), 2)*priorBatch.';

        pen = ((E + 1)./(O + 1)).^theta;
        R(:, idx) = base(:, idx).*(pen*phi(:, idx));
        R(:, idx) = R(:, idx)./(sum(R(:, idx), 1) + eps);

        O = O + R(:, idx)*phi(:, idx).';
        E = E + sum(R(:, idx), 2)*priorBatch.';
    end

    if norm(R - Rprev, "fro")/(norm(Rprev, "fro") + eps) < 1e-5
        break;
    end
end
end


function R = i_softassign(dist, sigma)
% Plain softmax over clusters, used to seed the responsibilities.
R = exp(-(dist - min(dist, [], 1))/sigma);
R = R./(sum(R, 1) + eps);
end


function sigma = i_autosigma(dist, target)
% Bisect on the bandwidth until a cell spreads over TARGET clusters on
% average. Perplexity (the exponential of the assignment entropy) rises
% monotonically with the bandwidth, so a plain bisection is enough.
lo = 1e-8;
hi = 1e8;
sigma = 1;
for iter = 1:50
    sigma = sqrt(lo*hi);
    if i_meanperplexity(dist, sigma) > target
        hi = sigma;
    else
        lo = sigma;
    end
end
end


function p = i_meanperplexity(dist, sigma)
logits = -dist/sigma;
logits = logits - max(logits, [], 1);
P = exp(logits);
P = P./(sum(P, 1) + eps);
p = mean(exp(-sum(P.*log(P + eps), 1)));
end


function Y = i_normalizecols(Y)
% L2-normalise each column, leaving all-zero columns alone.
nrm = vecnorm(Y, 2, 1);
nrm(nrm < eps) = 1;
Y = Y./nrm;
end
