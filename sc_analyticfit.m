function [T, xyzFit, curve, params] = sc_analyticfit(X, genelist, opts)
%SC_ANALYTICFIT  Closed-form replacement for the Spline-DV 3-D spline curve.
%
%   [T, xyzFit, curve, params] = SC_ANALYTICFIT(X, genelist) fits the
%   analytic mean / CV / dropout-rate curve implied by a gamma-Poisson
%   (negative-binomial) sampling model, instead of fitting a smoothing
%   spline through the gene cloud as SC_SPLINEFIT does.
%
%   With c the normalization scale factor (counts per c, default 1e4) and
%   L_j the raw library size of cell j, the curve is parameterized by the
%   normalized mean mu > 0:
%
%       x(mu) = log(1 + mu)
%       y(mu) = log(1 + sqrt(alpha/mu + phi))
%       z(mu) = (1/n) * sum_j (1 + phiZ * L_j * mu / c)^(-1/phiZ)
%
%   where
%       alpha = c * mean(1./L_j)   is fixed by the library sizes (no fitting)
%       phi                        is the gene-level overdispersion
%       phiZ -> 0 gives the Poisson limit  z(mu) = (1/n) sum_j exp(-L_j*mu/c),
%                                  i.e. the empirical Laplace transform of the
%                                  library sizes.
%
%   INPUTS
%     X          m-by-n matrix. Raw UMI counts by default; set
%                IsNormalized=true if X has already been library-size
%                normalized (LibSize is then required).
%     genelist   m-by-1 string of gene names (optional).
%
%   NAME-VALUE OPTIONS
%     LibSize         n-by-1 raw library sizes. Default sum(X,1).'
%     ScaleFactor     normalization scale c. Default 1e4, matching
%                     PKG.NORM_LIBSIZE.
%     Dispersion      [] (fit robustly, default), or a scalar phi; 0 = Poisson.
%     DropDispersion  overdispersion used for the dropout coordinate. Default
%                     0 (pure Poisson dropout); [] reuses Dispersion.
%     IsNormalized    Default false.
%     SortIt          sort T by deviation, descending. Default true.
%     OneSided        divide the deviation of genes below the curve by 100,
%                     as SC_SPLINEFIT does. Default true.
%     GridSize        resolution of the curve grid used to seed the
%                     projection. Default 20000. The projection is refined
%                     to sub-grid precision, so this only affects speed.
%     PlotIt          Default false.
%
%   OUTPUTS
%     T        table of genes, lgu, lgcv, dropr, d, pval, fdr, nearidx.
%              nearidx is 1:m so that xyzFit(T.nearidx,:) returns each gene's
%              foot point, preserving the call pattern used with
%              SC_SPLINEFIT.
%     xyzFit   m-by-3 foot point on the curve for each gene.
%     curve    struct of function handles x(mu), y(mu), z(mu), plus muHat,
%              the curve parameter at each gene's foot point.
%     params   struct with alpha, phi, phiZ, c, libSize.
%
%   See also SC_SPLINEFIT, SC_GENESTAT.

arguments
    X {mustBeNumeric, mustBeNonempty}
    genelist (:,1) string = strings(0,1)
    opts.LibSize (:,1) double = []
    opts.ScaleFactor (1,1) double {mustBePositive} = 1e4
    opts.Dispersion double {mustBeScalarOrEmpty, mustBeNonnegative} = []
    opts.DropDispersion double {mustBeScalarOrEmpty, mustBeNonnegative} = 0
    opts.IsNormalized (1,1) logical = false
    opts.SortIt (1,1) logical = true
    opts.OneSided (1,1) logical = true
    opts.GridSize (1,1) double {mustBePositive} = 20000
    opts.PlotIt (1,1) logical = false
end

% DOUBLE but not FULL: every step below preserves sparsity, and SCE.X is
% stored sparse single for all but the densest data. Densifying here cost
% 10 GB on a 25000 x 50000 matrix that is 1 GB sparse, which put HVG
% selection out of reach on exactly the datasets that need it.
X = double(X);
c = opts.ScaleFactor;

if isempty(genelist)
    genelist = "gene_" + string(1:size(X,1)).';
end

% ---- library sizes and normalization -----------------------------------
if isempty(opts.LibSize)
    if opts.IsNormalized
        error("sc_analyticfit:NeedLibSize", ...
            "LibSize is required when IsNormalized is true.");
    end
    libSize = sum(X, 1, "omitnan").';
else
    libSize = opts.LibSize;
end
if numel(libSize) ~= size(X, 2)
    error("sc_analyticfit:LibSizeLength", ...
        "LibSize must have one entry per cell.");
end

if opts.IsNormalized
    Xn = X;
else
    Xn = (X ./ libSize.') * c;
end

% ---- per-gene statistics (same definitions as SC_GENESTAT) -------------
keep = any(Xn > 0, 2);
if ~all(keep)
    Xn(~keep, :) = [];
    genelist(~keep) = [];
end
% FULL on the per-gene vectors, not on XN: these are m-long, so densifying
% them is cheap, and it keeps the table columns and the projection below
% dense whatever XN is. The values are identical either way.
% PKG.E_ROWVAR, not STD/VAR along dim 2: on a sparse matrix those walk
% every zero (7.0 s against 0.23 s at 20000 genes x 30000 cells); the
% values agree to ~1e-12 relative, and a constant row still gives 0.
u     = full(mean(Xn, 2, "omitnan"));
cv    = sqrt(pkg.e_rowvar(Xn, "omitnan")) ./ u;
lgu   = log1p(u);
lgcv  = log1p(cv);
dropr = 1 - full(sum(Xn > 0, 2)) ./ size(Xn, 2);

[~, order] = sortrows([lgu, dropr, lgcv], [1 3 2]);
lgu = lgu(order); lgcv = lgcv(order); dropr = dropr(order);
u = u(order); genes = genelist(order);

% ---- curve parameters --------------------------------------------------
alpha = c * mean(1 ./ libSize);          % closed form, nothing fitted

if isempty(opts.Dispersion)
    good = isfinite(u) & isfinite(lgcv) & u > 0;
    % robust (L1) fit of the single overdispersion parameter
    obj = @(lp) sum(abs(log1p(sqrt(alpha ./ u(good) + 10^lp)) - lgcv(good)));
    phi = 10 ^ fminbnd(obj, -6, 1);
else
    phi = opts.Dispersion;
end
phiZ = opts.DropDispersion;
if isempty(phiZ), phiZ = phi; end

curve.x = @(mu) log1p(mu);
curve.y = @(mu) log1p(sqrt(alpha ./ mu + phi));
curve.z = @(mu) i_dropout(mu, libSize, c, phiZ);

% ---- project every gene onto the curve ---------------------------------
P = [lgu, lgcv, dropr];
muLo  = max(min(u(u > 0)) / 10, realmin);
muHi  = max(u) * 10;
lgrid = linspace(log(muLo), log(muHi), opts.GridSize).';
[muHat, xyzFit] = i_project(P, lgrid, curve);

v = P - xyzFit;
d = vecnorm(v, 2, 2);
belowCurve = false(size(d));
if opts.OneSided
    belowCurve = v(:,2) < 0;
    d(belowCurve) = d(belowCurve) ./ 100;    % genes below the curve
end

% ---- significance, following SC_SPLINEFIT ------------------------------
% Deflating a distance by 100 is a heuristic "not a candidate" flag, not a
% transformation, so those genes are left out of the null and given p = 1
% rather than being allowed to define the scale everything else is judged
% against. Note there is no out-of-range class here: the analytic curve is
% defined for every mean, so only the below-curve genes are excluded, where
% SC_SPLINEFIT must also exclude genes past either end of its fit.
isCandidate = ~belowCurve;
pval = ones(size(d));
if any(isCandidate)
    pval(isCandidate) = pkg.e_deviationpvalue(d(isCandidate));
end
fdr = pkg.e_fdr(pval);

nearidx = (1:numel(d)).';
T = table(genes, lgu, lgcv, dropr, d, pval, fdr, nearidx);

curve.muHat = muHat;
params = struct("alpha", alpha, "phi", phi, "phiZ", phiZ, ...
    "c", c, "libSize", libSize);

if opts.SortIt
    T = sortrows(T, "d", "descend");
end

if opts.PlotIt
    mus = exp(linspace(log(min(u(u > 0))), log(max(u)), 500)).';
    figure;
    scatter3(P(:,1), P(:,2), P(:,3), "filled", "MarkerFaceAlpha", .1);
    hold on
    plot3(curve.x(mus), curve.y(mus), curve.z(mus), "-", "LineWidth", 4);
    xlabel("Mean, log"); ylabel("CV, log"); zlabel("Dropout rate");
    hold off
end
end

% =========================================================================
function z = i_dropout(mu, libSize, c, phi)
% Predicted fraction of zero counts at normalized mean MU: the average over
% cells of each cell's zero probability, (1 + phi*t)^(-1/phi) with
% t = libSize*mu/c (exp(-t) in the Poisson limit).
%
% It is evaluated at every grid point and again at every gene's projection
% -- about 40000 means times every cell -- and was most of this function's
% run time. Two exact rewrites, agreeing with the old form to ~1e-15:
%   - cells sharing a library size share a term, so the mean is a weighted
%     sum over the distinct sizes (one matrix-vector product);
%   - (1 + phi*t)^(-1/phi) is computed as exp(-log1p(phi*t)/phi), three
%     times faster than the non-integer power.
mu  = mu(:).';
z   = zeros(numel(mu), 1);
[s, ~, ic] = unique(libSize(:));
w   = accumarray(ic, 1).' / numel(libSize);    % share of cells at each size
blk = max(1, floor(2e7 / numel(s)));
for i = 1:blk:numel(mu)
    j = i:min(i + blk - 1, numel(mu));
    t = s * (mu(j) / c);
    if phi <= 1e-10
        z(j) = (w * exp(-t)).';
    else
        z(j) = (w * exp(-log1p(phi * t) / phi)).';
    end
end
end

% =========================================================================
function [muHat, F] = i_project(P, lgrid, curve)
% Nearest point on the curve: grid search, then one parabolic refinement in
% log(mu). The curve is smooth, so a single step reaches sub-grid precision.
mug = exp(lgrid);
C   = [curve.x(mug), curve.y(mug), curve.z(mug)];
% KNNSEARCH (a kd-tree) rather than DSEARCHN: the same nearest grid point
% for every gene, NaN rows included, 1.25 s against 0.07 s for 20000 genes
% on the default 20000-point grid.
k   = knnsearch(C, P);

edge = k <= 1 | k >= numel(lgrid);
k0   = min(max(k, 2), numel(lgrid) - 1);
l0 = lgrid(k0); lp = lgrid(k0 + 1);
dm = sum((C(k0-1,:) - P).^2, 2);
d0 = sum((C(k0,  :) - P).^2, 2);
dp = sum((C(k0+1,:) - P).^2, 2);

den  = dm - 2*d0 + dp;
step = zeros(size(den));
ok   = den > 0;
step(ok) = 0.5 * (dm(ok) - dp(ok)) ./ den(ok);   % in grid-spacing units
step = max(min(step, 1), -1);
lHat = l0 + step .* (lp - l0);
lHat(edge) = lgrid(k(edge));

muHat = exp(lHat);
F = [curve.x(muHat), curve.y(muHat), curve.z(muHat)];

worse = sum((F - P).^2, 2) > d0;                 % refinement must not hurt
if any(worse)
    muHat(worse)  = mug(k0(worse));
    F(worse, :)   = C(k0(worse), :);
end
end
