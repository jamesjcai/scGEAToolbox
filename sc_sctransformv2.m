function [R, params, Xcorrected] = sc_sctransformv2(X, opts)
% SC_SCTRANSFORMV2  Pure-MATLAB port of sctransform::vst(vst.flavor = "v2").
%
% Faithful reproduction of Seurat/sctransform's SCTransform v2 Pearson
% residuals: offset negative-binomial model with fixed slope, od_factor
% dispersion regularization by kernel regression against log10(gmean),
% Sheather-Jones bandwidth, min-variance floor and residual clipping.
%
% The only non-exact component is the per-gene step-1 fit: R uses the
% compiled glmGamPoi Gamma-Poisson estimator; here we use a MATLAB NB MLE.
%
% INPUT
%   X : G x C raw UMI counts (genes x cells)
%
% NAME-VALUE OPTIONS (sctransform::vst's step-1 subsampling)
%   NumGenes : genes whose NB model is fitted in step 1, drawn by density of
%              log10 geometric mean so every expression level is represented
%              (vst's n_genes, default 2000). The regression then carries the
%              fit to every gene.
%   NumCells : cells used for the step-1 fits (vst's n_cells, which
%              vst.flavor="v2" sets to 2000). Inf for both fits every
%              overdispersed gene on every cell, as this port did before.
%   Seed     : seed of the private random stream the two samples are drawn
%              from (default 1448145, Seurat's seed.use). The global stream
%              is not touched, and a given seed gives the same result.
%   On data with at most NumGenes step-1 genes and NumCells cells nothing is
%   sampled and the result is what fitting everything gives.
% OUTPUT
%   R          : G x C Pearson residuals (genes failing min_cells -> NaN row)
%   params     : struct with intermediate quantities (theta, intercept, mu, ...)
%   Xcorrected : G x C corrected counts, every cell re-expressed at the
%                median sequencing depth (sctransform::correct_counts, which
%                is Seurat's SCT "counts" slot; its "data" slot is log1p of
%                this). Computed only when a third output is requested.
%
% Checked against R sctransform 0.4.3 with vst.flavor="v2" on two simulated
% datasets: Pearson residuals agree to a per-gene correlation of 0.999996
% and 0.999999 (median), the top-200 genes by residual variance overlap
% 195/200, and corrected counts match to 99.7% of entries exactly. The
% Poisson/overdispersed call agrees on every gene. Theta is the loosest
% part, at a log10 correlation of 0.967, for the reason below.
% The step-1 subsampling (NumGenes, NumCells) came later, following vst.
% Against fitting every gene on every cell it gave, on the 7243 x 8260
% example, 196/200 of the top-200 genes by residual variance and a median
% per-gene residual correlation of 0.999985, in 6.2 s instead of 88 s.
%
% This replaced two earlier attempts, sctransform.m and sctransform_fast.m,
% which regularized only theta and did so with a smoothing spline rather
% than the kernel regression sctransform uses. Both looked close on the
% obvious measure -- per-gene correlation with R of 0.9999 and 0.9992 -- and
% both were useless for the job the transform is usually put to: of the top
% 200 genes R ranks by residual variance they recovered 2 and 3. Per-gene
% correlation is dominated by the counts and hides an error in the residual
% scale, which is exactly what regularizing the intercept as well as the
% dispersion fixes. Judge a change here on gene ranking, not correlation.

arguments
    X {mustBeNumeric}
    opts.NumGenes (1,1) double {mustBePositive} = 2000
    opts.NumCells (1,1) double {mustBePositive} = 2000
    opts.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1448145
end

X = double(full(X));
[G, C] = size(X);

min_cells = 5;
bw_adjust = 3;

% ----- cell / gene statistics --------------------------------------------
umi = sum(X, 1);                       % 1 x C total counts per cell
log_umi = log10(umi);                  % latent variable
mean_cell_sum = mean(umi);

amean = mean(X, 2);                    % G x 1 arithmetic mean
gvar  = var(X, 0, 2);                  % sample variance (N-1), matches rowVars
gmean = exp(mean(log(X + 1), 2)) - 1;  % row_gmean(eps = 1)
lgm   = log10(gmean);                  % log10 geometric mean
odfac = gvar - amean;                  % overdispersion factor (moment based)

% ----- gene filtering (min_cells) ----------------------------------------
genes_cell_count = sum(X >= 0.01, 2);
kept = genes_cell_count >= min_cells;  % genes carried through vst

% Poisson genes (get analytic offset, theta = Inf), among kept genes
low_mean = amean < 0.001;
all_poisson = kept & ((odfac <= 0) | low_mean);

% Step-1 training genes, as sctransform::vst picks them: a sample of cells,
% the genes detected in at least MIN_CELLS of those, a density-weighted
% sample of those genes, and finally only the overdispersed ones. The
% geometric mean and the overdispersion test use all cells, as vst does;
% only the model fits see the sampled cells. Fitting every overdispersed
% gene on every cell -- all this port did before -- took 96 s of 103 s on
% the 7243 x 8260 example.
stream = RandStream('twister', 'Seed', opts.Seed);
if opts.NumCells < C
    cells1 = sort(randperm(stream, C, opts.NumCells));
    genes1 = sum(X(:, cells1) > 0, 2) >= min_cells;
else
    cells1 = 1:C;
    genes1 = kept;
end
if opts.NumGenes < nnz(genes1)
    cand = find(genes1);
    pick = datasample(stream, cand, opts.NumGenes, 'Replace', false, ...
        'Weights', i_densityweights(lgm(cand)));
    genes1 = false(G, 1);
    genes1(pick) = true;
end
step1 = genes1 & kept & (odfac > 0);

% ----- step 1: per-gene offset NB fit (glmGamPoi substitute) --------------
theta1 = nan(G, 1);
b0_1   = nan(G, 1);
idx1 = find(step1);
umi1 = umi(cells1);
X1 = X(:, cells1);
mean_cell_sum1 = mean(umi1);
amean1 = mean(X1, 2);                  % starting values only; equal to AMEAN
gvar1 = var(X1, 0, 2);                 % and GVAR when no cells are sampled
for t = 1:numel(idx1)
    g = idx1(t);
    [b0_1(g), theta1(g)] = fit_nb_offset(X1(g, :), umi1, amean1(g), gvar1(g), mean_cell_sum1);
end

% method-of-moments Poisson override (diff_theta < 1e-3 -> theta = Inf)
predicted_theta = amean.^2 ./ (gvar - amean);
diff_theta = predicted_theta ./ theta1;
theta1(step1 & (diff_theta < 1e-3)) = Inf;
min_theta = 1e-7;
theta1(theta1 < min_theta) = min_theta;

% ----- dispersion_par (od_factor) on step-1 genes ------------------------
% dispersion_par = log10(1 + 10^lgm / theta)
disp_par = log10(1 + (10.^lgm) ./ theta1);   % 0 where theta = Inf
intercept1 = b0_1;

% ----- assemble training set for regularization --------------------------
% Remove outliers (per column) and theta=Inf and all_poisson genes.
train = step1;
lgm_s1 = lgm(train);
cols = [disp_par(train), intercept1(train), repmat(log(10), sum(train), 1)];
outlier = false(sum(train), 1);
for c = 1:size(cols, 2)
    outlier = outlier | is_outlier(cols(:, c), lgm_s1);
end
outlier = outlier | ~isfinite(theta1(train));      % exclude_poisson
tr_idx = find(train);
train(tr_idx(outlier)) = false;
% also drop all_poisson genes from training
train = train & ~all_poisson;

lgm_train = lgm(train);
disp_train = disp_par(train);
int_train  = intercept1(train);

% ----- kernel regression (ksmooth "normal") ------------------------------
bw = bw_sj(lgm_train) * bw_adjust;

kept_idx = find(kept);
lgm_kept = lgm(kept_idx);
xq = min(max(lgm_kept, min(lgm_train)), max(lgm_train));   % clamp to range

disp_fit = ksmooth_normal(lgm_train, disp_train, xq, bw);
int_fit  = ksmooth_normal(lgm_train, int_train,  xq, bw);

% ----- convert od_factor back to theta -----------------------------------
% for all_poisson genes dispersion_par is forced to 0 first
is_pois_kept = all_poisson(kept_idx);
disp_fit(is_pois_kept) = 0;
theta_fit = (10.^lgm_kept) ./ (10.^disp_fit - 1);   % Inf where disp_fit == 0

% Poisson genes: analytic offset intercept, theta = Inf
int_fit(is_pois_kept) = log(amean(kept_idx(is_pois_kept))) - log(mean_cell_sum);
theta_fit(is_pois_kept) = Inf;

% ----- residuals ----------------------------------------------------------
min_var = (median(X(X > 0)) / 5)^2;    % min_variance = "umi_median"
sn = sqrt(C);

R = nan(G, C);
Xk = X(kept_idx, :);
mu = exp(int_fit + log(10) .* log_umi);       % offset model, slope fixed
model_var = mu + mu.^2 ./ theta_fit;
model_var(model_var < min_var) = min_var;
Rk = (Xk - mu) ./ sqrt(model_var);
Rk(Rk >  sn) =  sn;
Rk(Rk < -sn) = -sn;
R(kept_idx, :) = Rk;

params = struct('theta_fit', theta_fit, 'intercept_fit', int_fit, ...
    'mu', mu, 'kept_idx', kept_idx, 'lgm', lgm, 'bw', bw, ...
    'min_var', min_var, 'theta1', theta1, 'intercept1', intercept1, ...
    'step1_genes', idx1, 'step1_cells', cells1);

% ----- corrected counts (sctransform::correct_counts) ---------------------
% Only on request: it is a second dense matrix the size of the input.
% Every cell is re-expressed at one common sequencing depth -- the median
% log UMI -- by reading the residual back through the model at that depth.
% That is what makes counts from cells sequenced to different depths
% comparable, and it is Seurat's SCT "counts" slot; its "data" slot is
% log1p of this. Reproduces R's correct_counts to 99.7% of entries exactly.
if nargout > 2
    Xcorrected = zeros(G, C);
    mu0 = exp(int_fit + log(10)*median(log_umi));
    var0 = mu0 + mu0.^2 ./ theta_fit;
    var0(var0 < min_var) = min_var;
    Xcorrected(kept_idx, :) = round(max(mu0 + Rk.*sqrt(var0), 0));
end
end

% =========================================================================
function [b0, theta] = fit_nb_offset(y, umi, am, gv, mean_cell_sum)
% Intercept-only NB with offset log(umi): mu = exp(b0) * umi.
% Estimate theta by maximizing the Cox-Reid adjusted profile likelihood
% (as glmGamPoi / edgeR / DESeq2 do); the plain MLE underestimates
% dispersion. Intercept b0 is the NB score-equation solution given theta.
y = y(:);
umi = umi(:);
b0_init = log(am) - log(mean_cell_sum);
if gv > am
    th_init = am^2 / (gv - am);
else
    th_init = 100;
end
th_init = min(max(th_init, 1e-3), 1e4);
opts = optimset('Display', 'off', 'TolX', 1e-6);
try
    lt = fminsearch(@(lt) -cr_apl(lt, y, umi, b0_init), log(th_init), opts);
    theta = exp(lt);
    b0 = solve_b0(theta, y, umi, b0_init);
catch
    b0 = b0_init;
    theta = th_init;
end
end

% -------------------------------------------------------------------------
function val = cr_apl(lt, y, umi, b0_init)
% Cox-Reid adjusted profile log-likelihood at log theta = lt.
theta = exp(lt);
b0 = solve_b0(theta, y, umi, b0_init);
mu = exp(b0) .* umi;
ll = sum(gammaln(y + theta) - gammaln(theta) - gammaln(y + 1) ...
    + theta .* (log(theta) - log(theta + mu)) ...
    + y .* (log(mu) - log(theta + mu)));
w = mu .* theta ./ (theta + mu);        % IRLS weight
val = ll - 0.5 * log(sum(w));           % intercept-only: |X'WX| = sum(w)
end

% -------------------------------------------------------------------------
function b0 = solve_b0(theta, y, umi, b0)
% Newton solve of the NB score equation for the intercept (offset log umi).
for it = 1:50
    mu = exp(b0) .* umi;
    S  = sum((y - mu) .* theta ./ (theta + mu));
    dS = sum(-mu .* theta ./ (theta + mu) ...
         - (y - mu) .* mu .* theta ./ (theta + mu).^2);
    step = S / dS;
    b0 = b0 - step;
    if abs(step) < 1e-10, break; end
end
end

% -------------------------------------------------------------------------
function yhat = ksmooth_normal(x, y, xout, bandwidth)
% Reproduce R's ksmooth(kernel = "normal"): the kernel is scaled so its
% quartiles are at +/- 0.25 * bandwidth, i.e. sigma = 0.3706506 * bandwidth.
sigma = 0.3706506 * bandwidth;
xout = xout(:);
x = x(:);
y = y(:);
yhat = zeros(numel(xout), 1);
for i = 1:numel(xout)
    w = exp(-0.5 * ((xout(i) - x) / sigma).^2);
    sw = sum(w);
    if sw > 0
        yhat(i) = sum(w .* y) / sw;
    else
        yhat(i) = NaN;
    end
end
end

% -------------------------------------------------------------------------
function out = is_outlier(y, x, th)
% Port of sctransform:::is_outlier.
if nargin < 3, th = 10; end
bin_width = (max(x) - min(x)) * bw_sj(x) / 2;
eps10 = eps * 10;
breaks1 = (min(x) - eps10):bin_width:(max(x) + bin_width);
breaks2 = (min(x) - eps10 - bin_width/2):bin_width:(max(x) + bin_width);
s1 = robust_scale_binned(y, x, breaks1);
s2 = robust_scale_binned(y, x, breaks2);
out = min(abs(s1), abs(s2)) > th;
end

% -------------------------------------------------------------------------
function score = robust_scale_binned(y, x, breaks)
bins = discretize(x, breaks);
score = zeros(numel(x), 1);
ub = unique(bins(~isnan(bins)));
for k = 1:numel(ub)
    m = bins == ub(k);
    score(m) = robust_scale(y(m));
end
end

% -------------------------------------------------------------------------
function s = robust_scale(x)
s = (x - median(x)) / (mad(x, 1) * 1.4826 + eps);
end

% -------------------------------------------------------------------------
function w = i_densityweights(x)
% sctransform's step-1 gene sampling weights, 1/(density of X + eps): R's
% density(x, bw = "nrd") on its default 512-point grid (cut = 3),
% interpolated back to each gene.
x = x(:);
bw = 1.06 * min(std(x), iqr7(x) / 1.34) * numel(x)^(-1/5);
grid = linspace(min(x) - 3*bw, max(x) + 3*bw, 512);
dens = ksdensity(x, grid, 'Bandwidth', bw);
w = 1 ./ (interp1(grid, dens, x) + eps);
end

% =========================================================================
function h = bw_sj(x)
% Port of R stats::bw.SJ(x, method = "ste").
x = x(:);
n = numel(x);
nb = 1000;
[d, cnt] = band_den_bin(x, nb);

scale = min(std(x), iqr7(x) / 1.349);
a = 1.24 * scale * n^(-1/7);
b = 1.23 * scale * n^(-1/9);
c1 = 1 / (2 * sqrt(pi) * n);
TD = -bw_phi6(n, d, cnt, b);
if ~isfinite(TD) || TD <= 0
    error('bw_sj:TD', 'sample is too sparse to find TD');
end
alph2 = 1.357 * (bw_phi4(n, d, cnt, a) / TD)^(1/7);
fSD = @(h) (c1 / bw_phi4(n, d, cnt, alph2 * h^(5/7)))^(1/5) - h;

hmax = 1.144 * scale * n^(-1/5);
lower = 0.1 * hmax;
upper = hmax;
tol = 0.1 * lower;
itry = 1;
while fSD(lower) * fSD(upper) > 0
    if itry > 99
        error('bw_sj:root', 'no solution in the specified range of bandwidths');
    end
    if mod(itry, 2)
        upper = upper * 1.2;
    else
        lower = lower / 1.2;
    end
    itry = itry + 1;
end
h = fzero(fSD, [lower, upper], optimset('TolX', tol));
end

% -------------------------------------------------------------------------
function [dd, cnt] = band_den_bin(x, nb)
% Reproduce R C_bw_den: binned pairwise distance counts.
xmin = min(x); xmax = max(x);
rang = (xmax - xmin) * 1.01;
dd = rang / nb;
cnt = zeros(nb, 1);
ii = fix(x / dd);            % (int) truncation toward zero, as in C
n = numel(x);
for i = 2:n
    j = 1:(i-1);
    iij = abs(ii(i) - ii(j)) + 1;   % +1 for 1-based indexing
    for k = 1:numel(iij)
        cnt(iij(k)) = cnt(iij(k)) + 1;
    end
end
end

% -------------------------------------------------------------------------
function u = bw_phi4(n, d, cnt, h)
DELMAX = 1000;
nbin = numel(cnt);
s = 0.0;
for i = 0:(nbin-1)
    delta = (i * d / h)^2;
    if delta >= DELMAX, break; end
    term = exp(-delta/2) * (delta^2 - 6*delta + 3);
    s = s + term * cnt(i+1);
end
s = 2*s + n*3;   % diagonal
u = s / (n * (n-1) * h^5 * sqrt(2*pi));
end

% -------------------------------------------------------------------------
function u = bw_phi6(n, d, cnt, h)
DELMAX = 1000;
nbin = numel(cnt);
s = 0.0;
for i = 0:(nbin-1)
    delta = (i * d / h)^2;
    if delta >= DELMAX, break; end
    term = exp(-delta/2) * (delta^3 - 15*delta^2 + 45*delta - 15);
    s = s + term * cnt(i+1);
end
s = 2*s - n*15;
u = s / (n * (n-1) * h^7 * sqrt(2*pi));
end

% -------------------------------------------------------------------------
function q = iqr7(x)
% R-type-7 IQR to match R's IQR().
q = quantile7(x, 0.75) - quantile7(x, 0.25);
end

function v = quantile7(x, p)
x = sort(x(:));
n = numel(x);
h = (n - 1) * p + 1;
lo = floor(h);
hi = min(lo + 1, n);
v = x(lo) + (h - lo) * (x(hi) - x(lo));
end
