function [T, info] = sc_powerseek(X, genelist, options)
%SC_POWERSEEK Rank genes by how much they drive cell-state covariance structure.
%
%   T = sc_powerseek(X, genelist)
%   [T, info] = sc_powerseek(X, genelist, Name=Value)
%
%   Power-Seek (Ghosh, Chakrabarti & Raju, Cell Systems 2026,
%   doi:10.1016/j.cels.2026.101707). Cells that share a slowly relaxing
%   expression state are correlated with each other, which turns the
%   eigenvalue-vs-rank curve of the cell-cell covariance matrix into a power
%   law; independent cells give a flat Marchenko-Pastur spectrum instead.
%   Each gene is scored by how much the top (power-law) eigenvalues drop
%   when that gene alone is left out:
%
%       delta_g = sum_{j<=n} ( lambda_j - lambda~_j(g) )
%
%   with lambda the eigenvalues of C = Z*Z'/ng and lambda~(g) those of the
%   same matrix built without gene g, divided by ng-1. A gene with
%   delta_g > 0 holds up the cell-cell structure; most genes give a small
%   negative delta_g.
%
%   What the score measures is covariation with the dominant structure
%   among cells. The authors read high-delta genes as "memory" genes whose
%   expression is inherited across divisions, but a single snapshot cannot
%   tell lineage memory apart from any other cause of shared cell states,
%   such as a persistent transcriptional program or a common
%   microenvironment. Run it within one cell type: differences between cell
%   types dominate the spectrum otherwise.
%
%   Inputs
%     X         genes x cells raw counts (full or sparse)
%     genelist  gene names, one per row of X
%
%   Name-Value arguments
%     NumGenes        genes kept by coefficient of variation (default 3000;
%                     the paper uses 3000-5000)
%     MinCells        a gene must be detected in more than this many cells
%                     (default min(200, 5% of cells); the paper uses 200)
%     RemoveOutliers  drop genes whose log nonzero expression has IQR
%                     outliers, i.e. very high in a few cells (default true)
%     NumEigenvalues  size n of the power-law regime; NaN (default) picks the
%                     first peak of R^2 of the log-log fit against n, as
%                     the authors do, or the smallest fit (20) when R^2
%                     falls from there on
%     Preprocess      false skips library-size scaling, the gene filters
%                     and the CV selection, for X that is already normalized
%                     expression; genes are still z-scored (default true)
%     Shuffle         also compute the spectrum after permuting each gene
%                     across cells, the null without shared states
%                     (default true)
%     Seed            seed of the local stream used for the shuffle
%     ProgressFcn     handle called as f(fraction) while genes are scored
%
%   Outputs
%     T     table sorted by Delta, descending: Gene, Delta, Rank
%     info  struct: Eigenvalues, EigenvaluesShuffled, NumEigenvalues,
%           NumEigenvaluesTried, RSquared, HasPowerLaw, Exponent,
%           ExponentShuffled, NumCells, NumGenes. HasPowerLaw is false
%           when the log-log fit over the chosen n has R^2 below 0.85, or,
%           with Shuffle on, when its slope is under 1.5 times the shuffled
%           one: the top of the spectrum then shows no shared structure and
%           the ranking should not be trusted. Both cutoffs are this
%           implementation's, not the paper's. It is NaN when n was given.
%
%   The leave-one-gene-out spectra are exact. Removing gene g subtracts the
%   rank-one term x_g*x_g' from Z*Z', so its eigenvalues are the roots of a
%   secular equation built from one eigendecomposition, not ng new ones.
%
%   See also sc_hvg, sc_genestats

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    options.NumGenes (1, 1) double {mustBeInteger, mustBePositive} = 3000
    options.MinCells (1, 1) double = NaN
    options.RemoveOutliers (1, 1) logical = true
    options.NumEigenvalues (1, 1) double = NaN
    options.Preprocess (1, 1) logical = true
    options.Shuffle (1, 1) logical = true
    options.Seed (1, 1) double = 1
    options.ProgressFcn = []
end

if size(X, 1) ~= numel(genelist)
    error("sc_powerseek:sizeMismatch", ...
        "X has %d rows but genelist has %d names. Pass X as genes x cells.", ...
        size(X, 1), numel(genelist));
end

if options.Preprocess
    [Z, genes] = preprocessCounts(X, genelist, options);
else
    genes = genelist;
    Z = full(double(X))';
    Z = (Z - mean(Z, 1))./std(Z, 1, 1);
end
[nc, ng] = size(Z);
if ng < 3 || nc < 3
    error("sc_powerseek:tooFewGenes", ...
        "Only %d genes and %d cells pass the filters. Lower MinCells or use more cells.", ...
        ng, nc);
end

% Eigenvalues d of Z*Z' (unscaled) and W, whose row g holds the components
% of gene g's column in the eigenbasis: W(g,i)^2 = s_i^2 * V(g,i)^2. Both
% come from whichever Gram matrix is smaller.
[d, W] = gramEigen(Z);
ev = d/ng;

if isnan(options.NumEigenvalues)
    [n, nTried, r2, hasPowerLaw] = findPowerLawRegime(ev);
else
    n = options.NumEigenvalues;
    nTried = [];
    r2 = [];
    hasPowerLaw = NaN;
end
n = min(n, numel(d) - 1);
if n < 1
    error("sc_powerseek:noSpectrum", "The covariance matrix has too few nonzero eigenvalues.");
end

delta = leaveOneOutDelta(d, W.^2, n, ng, options.ProgressFcn);

[~, order] = sort(delta, "descend");
rank = zeros(ng, 1);
rank(order) = (1:ng)';
T = table(genes(order), delta(order), rank(order), ...
    VariableNames=["Gene", "Delta", "Rank"]);

info = struct();
info.Eigenvalues = ev;
info.EigenvaluesShuffled = [];
info.NumEigenvalues = n;
info.NumEigenvaluesTried = nTried;
info.HasPowerLaw = hasPowerLaw;
info.RSquared = r2;
info.Exponent = powerLawSlope(ev, n);
info.ExponentShuffled = NaN;
info.NumCells = nc;
info.NumGenes = ng;
if options.Shuffle
    stream = RandStream("mt19937ar", Seed=options.Seed);
    Zs = Z;
    for k = 1:ng
        Zs(:, k) = Z(randperm(stream, nc), k);
    end
    info.EigenvaluesShuffled = gramEigen(Zs)/ng;
    info.ExponentShuffled = powerLawSlope(info.EigenvaluesShuffled, n);
    % The top of a structureless (Marchenko-Pastur) spectrum also fits a
    % line on log-log axes; what sets shared structure apart is a slope
    % well beyond the shuffled one.
    minSlopeRatio = 1.5;
    if ~isnan(info.HasPowerLaw)
        info.HasPowerLaw = info.HasPowerLaw && ...
            info.Exponent <= minSlopeRatio*info.ExponentShuffled;
    end
end
end

function [Z, genes] = preprocessCounts(X, genelist, options)
% Library-size scaling without a log (a log transform weakens the
% correlations the method looks for), sparsity and outlier filters, top
% genes by coefficient of variation, then z-scores per gene.
X = double(X);
nCells = size(X, 2);
keep = full(sum(X ~= 0, 2)) > 0;
X = X(keep, :);
genelist = genelist(keep);

total = full(sum(X, 1));
hasCounts = total > 0;
X = X(:, hasCounts)./total(hasCounts);
nCells = min(nCells, nnz(hasCounts));

minCells = options.MinCells;
if isnan(minCells)
    minCells = min(200, ceil(0.05*nCells));
end
keep = full(sum(X ~= 0, 2)) > minCells;
X = X(keep, :);
genelist = genelist(keep);

if options.RemoveOutliers
    keep = ~hasExpressionOutliers(X);
    X = X(keep, :);
    genelist = genelist(keep);
end

mu = full(mean(X, 2));
sd = full(std(X, 0, 2));
cv = sd./mu;
[~, order] = sort(cv, "descend");
order = order(1:min(options.NumGenes, numel(order)));
genes = genelist(order);

Z = full(X(order, :))';
Z = (Z - mean(Z, 1))./std(Z, 1, 1);
end

function isOutlier = hasExpressionOutliers(X)
% A gene is flagged when any log nonzero value lies outside 1.5 x the
% 9th-91st percentile range beyond those percentiles.
lowerPercentile = 9;
upperPercentile = 91;
whisker = 1.5;
Xt = X';
nGenes = size(Xt, 2);
isOutlier = false(nGenes, 1);
for k = 1:nGenes
    v = nonzeros(Xt(:, k));
    v = log(full(v));
    q = prctile(v, [lowerPercentile, upperPercentile], Method="inclusive");
    spread = q(2) - q(1);
    isOutlier(k) = any((v < q(1) - whisker*spread) | (v > q(2) + whisker*spread));
end
end

function [d, W] = gramEigen(Z)
% Nonzero eigenvalues of Z*Z' in descending order; with nargout 2, also
% W = Z'*U, the gene loadings scaled by the singular values. Centering
% costs one rank, so at most min(nc - 1, ng) eigenvalues are kept.
[nc, ng] = size(Z);
if nc <= ng
    G = Z*Z';
else
    G = Z'*Z;
end
G = (G + G')/2;
maxRank = min(nc - 1, ng);
if nargout < 2
    d = sort(eig(G), "descend");
    d = d(1:maxRank);
    d = d(d > d(1)*1e-10);
    return
end
[V, D] = eig(G, "vector");
[d, order] = sort(D, "descend");
keep = order(1:maxRank);
d = d(1:maxRank);
keep = keep(d > d(1)*1e-10);
d = D(keep);
V = V(:, keep);
if nc <= ng
    W = Z'*V;
else
    W = V.*sqrt(d');
end
end

function [n, nTried, r2, hasPowerLaw] = findPowerLawRegime(ev)
% R^2 of a straight-line fit of log(eigenvalue) on log(rank) over the top
% m eigenvalues, for m = 20, 22, ...; the first local peak marks where the
% spectrum leaves the power law. R^2 rises again further out, as the noisy
% bulk fills in. When R^2 falls from the first fit on, the regime is
% shorter than that fit and the first peak is the noise peak, so the
% smallest fit is used instead. Either way the regime counts as a power law
% only if its fit reaches minRSquared (the paper's fits on real data had
% 0.89 and 0.92; the cutoff is this implementation's, not the paper's).
firstSize = 20;
stepSize = 2;
dipTolerance = 0.01;
minRSquared = 0.85;
nTried = (firstSize:stepSize:numel(ev) - 1)';
hasPowerLaw = false;
if numel(nTried) < 3
    n = max(1, min(firstSize, numel(ev) - 1));
    r2 = [];
    return
end
x = log((1:numel(ev))');
y = log(ev);
r2 = zeros(size(nTried));
for k = 1:numel(nTried)
    m = nTried(k);
    c = corrcoef(x(1:m), y(1:m));
    r2(k) = c(1, 2)^2;
end
slope = diff(r2);
peak = find((slope(1:end - 1) > 0) & (slope(2:end) < 0), 1) + 1;
% A dip of more than dipTolerance below the first fit before the peak
% means R^2 fell first.
if isempty(peak) || min(r2(1:peak)) < r2(1) - dipTolerance
    peak = 1;
end
n = nTried(peak);
hasPowerLaw = r2(peak) >= minRSquared;
end

function b = powerLawSlope(ev, n)
n = min(n, numel(ev));
p = polyfit(log((1:n)'), log(ev(1:n)), 1);
b = p(1);
end

function delta = leaveOneOutDelta(d, z2, n, ng, progressFcn)
% For each gene g, the top n eigenvalues mu of Z*Z' - x_g*x_g' solve
%   1 = sum_i z2(g,i)/(d_i - mu),
% one root in each interval (d_{j+1}, d_j). Solved for t = d_j - mu, so the
% small distance to the nearest pole keeps full precision. Each step is
% Gragg's fixed-weight step: the poles at either end of the interval are
% modelled exactly and the rest of each side by a matched constant plus a
% pole term, which turns the update into a quadratic with one root in the
% interval. A bracket guards it, stepping by secant or bisection when the
% model step leaves it.
maxIterations = 50;
tolerance = 1e-12;
maxElements = 8e6;                             % per genes x r x n work array

r = numel(d);
dTop = d(1:n)';                                % 1 x n
gap = d(1:n)' - d(2:n + 1)';                   % 1 x n, bracket width
D = reshape(d - dTop, [1, r, n]);              % d_i - d_j
isUpper = reshape((1:r)' <= (1:n), [1, r, n]); % poles at or left of t = 0

nGenes = size(z2, 1);
blockSize = max(1, floor(maxElements/(r*n)));
delta = zeros(nGenes, 1);
for first = 1:blockSize:nGenes
    rows = first:min(first + blockSize - 1, nGenes);
    nb = numel(rows);
    zz = z2(rows, :);                          % nb x r
    g = repmat(gap, nb, 1);
    lo = zeros(nb, n);
    hi = g;
    hLo = -inf(nb, n);                         % h at the bracket ends: the
    hHi = inf(nb, n);                          % poles, until replaced
    t = g/2;

    for iteration = 1:maxIterations
        q = 1./(D + reshape(t, [nb, 1, n]));
        terms = zz.*q;
        slopes = terms.*q;
        psiUpper = reshape(sum(terms.*isUpper, 2), [nb, n]);
        psiLower = reshape(sum(terms, 2), [nb, n]) - psiUpper;
        slopeUpper = reshape(sum(slopes.*isUpper, 2), [nb, n]);
        slopeLower = reshape(sum(slopes, 2), [nb, n]) - slopeUpper;

        h = 1 - psiUpper - psiLower;
        below = h < 0;
        above = h > 0;
        lo(below) = t(below);
        hLo(below) = h(below);
        hi(above) = t(above);
        hHi(above) = h(above);

        bUpper = slopeUpper.*t.^2;
        bLower = slopeLower.*(t - g).^2;
        c = h + bUpper./t + bLower./(t - g);
        B = c.*g + bUpper + bLower;
        tNew = 2*bUpper.*g./(B + sqrt(B.^2 - 4*c.*bUpper.*g));

        % A step that leaves the bracket falls back to the secant between
        % its ends, or to the midpoint while an end is still a pole.
        outside = ~(tNew >= lo & tNew <= hi);
        secant = lo - hLo.*(hi - lo)./(hHi - hLo);
        useSecant = outside & isfinite(hLo) & isfinite(hHi);
        tNew(useSecant) = secant(useSecant);
        bisect = outside & ~useSecant;
        tNew(bisect) = (lo(bisect) + hi(bisect))/2;

        % Judged on the residual, not the step: a secant step can creep.
        converged = abs(h) <= tolerance*(1 + psiUpper - psiLower) | ...
            hi - lo <= tolerance*t;
        t(~converged) = tNew(~converged);
        if all(converged, "all")
            break
        end
    end

    % lambda_j - lambda~_j = d_j/ng - (d_j - t_j)/(ng - 1)
    delta(rows) = sum(t/(ng - 1) - dTop/(ng*(ng - 1)), 2);
    if ~isempty(progressFcn)
        progressFcn(rows(end)/nGenes);
    end
end
end
