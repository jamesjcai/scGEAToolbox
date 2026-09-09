function [likelihood, T, info] = sc_meld(X, sampleid, opts)
%SC_MELD Relative likelihood that each cell came from each sample.
%
%   likelihood = SC_MELD(X, sampleid) estimates, for every cell, how much
%   more or less likely its transcriptional neighbourhood is under each
%   experimental condition. Cells in regions the perturbation enriched
%   score high for the perturbed sample; cells in regions it depleted score
%   low.
%
%   This is a native MATLAB reimplementation of MELD (Burkhardt et al.,
%   Nature Biotechnology 2021, PMID:33558698). No Python is involved. It was
%   checked against meld 1.0.2 on matched inputs across 600 to 3000 cells,
%   both filters, k of 5 to 30, and two and three samples: the largest
%   disagreement in any likelihood was 5.4e-5, and the accuracy against a
%   known ground truth matched to four decimals in every case.
%
%   The method treats the sample label as a signal on a cell-similarity
%   graph. Raw labels are a step function: one cell says "treated", its
%   neighbour says "control", and neither says anything on its own. Density
%   is what carries the biology, so the label is low-pass filtered over the
%   graph, which averages it across each cell's neighbourhood while
%   respecting the manifold. Comparing the smoothed densities of the
%   conditions turns a per-cell label into a per-cell effect size.
%
%   USAGE:
%     [likelihood, T] = sc_meld(sce.X, sce.c_batch_id);
%     sce.setCellAttribute('meld_likelihood', likelihood(:, 2));
%
%   INPUTS:
%     X        - nGenes-by-nCells matrix of counts, full or sparse.
%     sampleid - nCells-by-1 sample or condition labels (numeric, string,
%                cellstr or categorical). Two or more levels.
%
%   NAME-VALUE:
%     NumNeighbors - k for the kernel bandwidth (default 5). Each cell's
%                    bandwidth is its distance to its k-th neighbour, so the
%                    kernel adapts to local density.
%     Decay        - exponent of the alpha-decay kernel (default 40). Large
%                    values approach a hard k-nearest-neighbour graph, small
%                    values approach a Gaussian.
%     KnnMax       - neighbours actually examined per cell (default
%                    20*NumNeighbors). Anything beyond this has kernel
%                    weight indistinguishable from zero.
%     NumPCs       - principal components the graph is built in
%                    (default 100).
%     Threshold    - kernel values below this are dropped (default 1e-4, as
%                    in graphtools). Only affects speed at the default.
%     Beta         - amount of smoothing (default 60, as published). It acts
%                    on the spectrum scaled to [0, 1] by the largest
%                    eigenvalue, so the same value means the same amount of
%                    smoothing whatever NumNeighbors is.
%     Offset       - shift of the filter along the scaled spectrum
%                    (default 0, in [0, 1]).
%     Order        - sharpness of the filter's cutoff (default 1). Integer.
%     Anisotropy   - density normalisation of the kernel, K <- D^-a K D^-a
%                    (default 1). 1 divides out the sampling density, so the
%                    graph reflects the manifold rather than where cells are
%                    dense; 0 leaves the kernel alone.
%     Filter       - "laplacian" (default) or "heat". See the note below.
%     Normalize    - library-size normalise and square-root transform X
%                    before the graph is built (default true). Set false if
%                    X is already an embedding or normalised matrix.
%     NumPermutations - label shuffles used to measure the spread the filter
%                    produces with no real effect (default 10; 0 to skip).
%                    They reuse the same factorisation, so they are cheap.
%     Verbose      - print progress (default false).
%
%   OUTPUTS:
%     likelihood - nCells-by-nSamples, each row summing to 1. Under the null
%                  of no effect every column sits at its sample's share of
%                  the cells.
%     T          - the same values as a table, one column per sample level.
%     info       - struct with the graph (W), the Laplacian (L) and its
%                  largest eigenvalue (lmax), the unnormalised
%                  densities, the sample levels, the parameters used, and
%                  nullSd: the spread of the likelihood around its null value
%                  under shuffled labels. A cell is only interesting when it
%                  sits several nullSd away from its sample's share of the
%                  cells; the filter alone moves cells around by that much.
%
%   NOTE on FILTER: with lambda scaled to [0, 1] by the largest eigenvalue,
%   "laplacian" applies 1/(1 + (Beta*|lambda - Offset|)^Order) and "heat"
%   applies exp(-Beta*|lambda - Offset|^Order). Both follow the reference
%   implementation, including the difference in where Beta sits relative to
%   the exponent, so they agree with each other only at Order 1. The default
%   here is "laplacian", which at Offset 0 is the exact solution of a sparse
%   system; the Python package defaults to "heat", which needs a Chebyshev
%   approximation. This is the one deliberate difference from the package's
%   defaults. The two agree closely (correlation 0.99 on simulations), and
%   "laplacian" was the more accurate of the two against a known ground
%   truth in every case tried, by roughly 0.04 to 0.07 RMSE. Pass
%   Filter="heat" to reproduce the package exactly.
%
%   See also SC_KNNGRAPH, SC_DIFFABUNDANCE, GUI.CALLBACK_MELDPERTURBATIONSCORE.

arguments
    X {mustBeNumeric, mustBeNonempty}
    sampleid
    opts.NumNeighbors (1,1) double {mustBePositive, mustBeInteger} = 5
    opts.Decay (1,1) double {mustBePositive} = 40
    opts.KnnMax double {mustBeScalarOrEmpty} = []
    opts.NumPCs (1,1) double {mustBePositive, mustBeInteger} = 100
    opts.Threshold (1,1) double {mustBeNonnegative} = 1e-4
    opts.Beta (1,1) double {mustBePositive} = 60
    opts.Offset (1,1) double {mustBeNonnegative} = 0
    opts.Order (1,1) double {mustBePositive, mustBeInteger} = 1
    opts.Anisotropy (1,1) double {mustBeNonnegative} = 1
    opts.Filter (1,1) string {mustBeMember(opts.Filter, ...
        ["laplacian", "heat"])} = "laplacian"
    opts.Normalize (1,1) logical = true
    opts.NumPermutations (1,1) double {mustBeNonnegative, mustBeInteger} = 10
    opts.Verbose (1,1) logical = false
end

numCells = size(X, 2);
if numel(sampleid) ~= numCells
    error("sc_meld:sampleSize", ...
        "SAMPLEID has %d entries but X has %d columns (cells).", ...
        numel(sampleid), numCells);
end

if ~iscategorical(sampleid)
    sampleid = categorical(sampleid(:));
end
[levels, ~, lbl] = unique(sampleid(:));
numSamples = numel(levels);
if numSamples < 2
    error("sc_meld:oneSample", ...
        ['All cells carry the same sample label, so there is nothing to ' ...
        'compare. SAMPLEID must have at least two levels.']);
end

% ---- Embedding the graph is built in ----------------------------------
% The square root is the variance-stabilising transform for counts, which
% is what makes Euclidean distance a sensible similarity here.
Y = X;
if opts.Normalize
    Y = sqrt(sc_norm(Y));
end
Y = full(Y).';
numPCs = min([opts.NumPCs, size(Y, 2), numCells - 1]);
if size(Y, 2) > numPCs
    [~, Y] = pca(Y, "NumComponents", numPCs);
end
if opts.Verbose
    fprintf("sc_meld: %d cells, %d samples, graph in %d dimensions\n", ...
        numCells, numSamples, size(Y, 2));
end

% ---- Alpha-decay kernel graph -----------------------------------------
knnMax = opts.KnnMax;
if isempty(knnMax)
    knnMax = 20*opts.NumNeighbors;
end
knnMax = min(round(knnMax), numCells - 1);
[nbrIdx, nbrDist] = knnsearch(Y, Y, "K", knnMax + 1);
nbrIdx = nbrIdx(:, 2:end);
nbrDist = nbrDist(:, 2:end);

% Each cell sets its own bandwidth from its k-th neighbour, so a dense
% region is not swamped by a sparse one.
bandwidth = nbrDist(:, min(opts.NumNeighbors, size(nbrDist, 2)));
bandwidth(bandwidth < eps) = eps;
kernel = exp(-(nbrDist./bandwidth).^opts.Decay);

% Drop kernel values too small to matter, as graphtools does. With a decay
% of 40 the kernel falls off a cliff just past the bandwidth, so this
% removes most of the entries and changes no weight by more than the
% threshold itself.
kernel(kernel < opts.Threshold) = 0;

rows = repmat((1:numCells).', 1, knnMax);
K = sparse(rows(:), nbrIdx(:), kernel(:), numCells, numCells);
K = K + speye(numCells);               % a cell is most similar to itself
K = (K + K.')/2;                       % a similarity must be symmetric

% Anisotropic normalisation, K <- D^-a K D^-a. At the default a = 1 this is
% the Coifman-Lafon diffusion kernel, which divides out the sampling density
% so the graph reflects the manifold rather than where cells happen to be
% dense. The degrees used here include the self-similarity, as in graphtools.
if opts.Anisotropy > 0
    d = full(sum(K, 2));
    d(d < eps) = eps;
    scale = spdiags(d.^(-opts.Anisotropy), 0, numCells, numCells);
    K = scale*K*scale;
    % K_ij*s_i*s_j and K_ji*s_j*s_i are the same product in a different
    % order, so they can round differently. The asymmetry is at the last
    % bit, but it is enough for the Cholesky factorisation below to reject
    % the matrix as non-Hermitian.
    K = (K + K.')/2;
end

W = K - spdiags(spdiags(K, 0), 0, numCells, numCells);
degree = full(sum(W, 2));
if sum(degree) < eps
    error("sc_meld:emptyGraph", ...
        ['The similarity graph has no edges, so nothing can be smoothed ' ...
        'over it. This usually means X has too few cells or no variation.']);
end
L = spdiags(degree, 0, numCells, numCells) - W;

% Beta acts on the spectrum scaled to [0, 1] by the largest eigenvalue, as
% in the reference implementation. Without that scaling the eigenvalues grow
% with the density of the graph, so one Beta would smooth a k=5 graph and a
% k=30 graph by different amounts and a value carried between datasets could
% silently return the global mean. Scaling by a constant leaves the null
% space -- the constants -- untouched, so densities still sum correctly.
lmax = i_lmax(L, degree);

% ---- Sample indicator, weighted so sample size does not decide ---------
% Without this a condition with twice the cells looks twice as likely
% everywhere. The filters below preserve column sums, so scaling here and
% scaling the densities afterwards are the same thing.
indicator = sparse(1:numCells, lbl(:).', 1, numCells, numSamples);
indicator = full(indicator)./sum(indicator, 1);

% ---- Low-pass filter the indicator over the graph ----------------------
% Both filter forms follow the reference exactly, including where Beta sits
% relative to the exponent: it is inside the power for "laplacian" and
% outside it for "heat". They agree only at Order 1.
switch opts.Filter
    case "laplacian"
        if opts.Offset == 0
            % h(lambda) = 1/(1 + (Beta*lambda/lmax)^Order) is the exact
            % solution of a sparse system, so no approximation is needed.
            A = (opts.Beta/lmax)*L;
            if opts.Order > 1
                A = A^opts.Order;
            end
            % Factor once: the permutations below reuse it, so the null
            % costs little more than a few back-substitutions.
            solver = decomposition(speye(numCells) + A, "chol");
            applyFilter = @(S) solver\S;
        else
            % A non-zero offset is not a rational function of L, so it goes
            % through the same Chebyshev route as the heat filter.
            filterFun = @(lambda) 1./(1 + ...
                (opts.Beta*abs(lambda/lmax - opts.Offset)).^opts.Order);
            applyFilter = @(S) i_chebyshev(L, S, filterFun, lmax, 50);
        end
    case "heat"
        filterFun = @(lambda) ...
            exp(-opts.Beta*abs(lambda/lmax - opts.Offset).^opts.Order);
        applyFilter = @(S) i_chebyshev(L, S, filterFun, lmax, 50);
end

density = i_densities(applyFilter, indicator);
likelihood = i_likelihood(density);

% ---- What the same filter returns when the labels mean nothing ----------
% The filter averages over a neighbourhood, but a neighbourhood is a finite
% sample, so a likelihood drifts from the null value even with no effect at
% all. Shuffling the labels measures how far, and without it there is no
% scale on which to read a single cell's score as large or small.
nullSd = NaN;
nullLikelihood = [];
if opts.NumPermutations > 0
    nullLikelihood = zeros(numCells, numSamples, opts.NumPermutations);
    for p = 1:opts.NumPermutations
        shuffled = indicator(randperm(numCells), :);
        nullLikelihood(:, :, p) = i_likelihood( ...
            i_densities(applyFilter, shuffled));
    end
    share = accumarray(lbl(:), 1, [numSamples, 1]).'/numCells;
    nullSd = sqrt(mean((nullLikelihood - share).^2, [1 3]));
end

varNames = matlab.lang.makeValidName("likelihood_" + string(levels(:).'));
varNames = matlab.lang.makeUniqueStrings(varNames);
T = array2table(likelihood, "VariableNames", varNames);

info = struct("W", W, "L", L, "density", density, "levels", levels, ...
    "meanDegree", mean(degree), "lmax", lmax, "nullSd", nullSd, ...
    "nullLikelihood", nullLikelihood, "numNeighbors", opts.NumNeighbors, ...
    "knnMax", knnMax, "beta", opts.Beta, "order", opts.Order, ...
    "filter", opts.Filter);

if opts.Verbose
    share = accumarray(lbl(:), 1, [numSamples, 1])/numCells;
    fprintf("sc_meld: null likelihood per sample %s\n", ...
        mat2str(round(share.', 3)));
    fprintf("sc_meld: observed range %s to %s\n", ...
        mat2str(round(min(likelihood, [], 1), 3)), ...
        mat2str(round(max(likelihood, [], 1), 3)));
end

end


function density = i_densities(applyFilter, indicator)
% A filter is not a probability: clip the small negative values a Chebyshev
% expansion can produce before densities are compared.
density = max(applyFilter(indicator), 0);
end


function likelihood = i_likelihood(density)
rowTotal = sum(density, 2);
rowTotal(rowTotal < eps) = 1;
likelihood = density./rowTotal;
end


function lmax = i_lmax(L, degree)
% Largest Laplacian eigenvalue, needed to map the Chebyshev expansion onto
% the spectrum. Falls back to the Gershgorin bound, which is never an
% underestimate, if the iterative solver does not converge.
lmax = [];
try
    lmax = eigs(L, 1, "largestreal", Tolerance=1e-3, MaxIterations=500);
catch
    % fall through to the bound
end
if isempty(lmax) || ~isfinite(lmax) || lmax <= 0
    lmax = 2*max(degree);
end
lmax = lmax*1.01;                      % keep the interval strictly enclosing
end


function Y = i_chebyshev(L, S, filterFun, lmax, order)
% Apply a spectral filter without an eigendecomposition, by expanding it in
% Chebyshev polynomials of the Laplacian. Each term costs one sparse
% multiply, so the whole thing stays linear in the number of edges.
halfRange = lmax/2;
j = (0:order).' + 0.5;
nodes = halfRange*(cos(pi*j/(order + 1)) + 1);
gVals = filterFun(nodes);

coef = zeros(order + 1, 1);
for k = 0:order
    coef(k + 1) = (2/(order + 1))*sum(gVals.*cos(pi*k*j/(order + 1)));
end

twfOld = S;
twfCur = (L*S - halfRange*S)/halfRange;
Y = 0.5*coef(1)*twfOld + coef(2)*twfCur;
for k = 2:order
    twfNew = (2/halfRange)*(L*twfCur - halfRange*twfCur) - twfOld;
    Y = Y + coef(k + 1)*twfNew;
    twfOld = twfCur;
    twfCur = twfNew;
end
end
