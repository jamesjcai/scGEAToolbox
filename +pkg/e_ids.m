function [C, p] = e_ids(X, Y, options)
%E_IDS InterDependence Score between the columns of X (and optionally Y).
%
%   C = pkg.e_ids(X)
%   C = pkg.e_ids(X, Y)
%   [C, p] = pkg.e_ids(___, NumPermutations=100)
%   C = pkg.e_ids(___, PNorm="max", NumTerms=6, Bandwidth=0.5)
%
%   X is n-by-dx with samples in rows and variables in columns, the same
%   orientation as CORR. With X alone, C is dx-by-dx; with Y (n-by-dy), C is
%   dx-by-dy. Values lie in [0, 1]; 0 indicates independence.
%
%   Each variable v is lifted into NumTerms features
%       phi_i(v) = exp(-Bandwidth*v.^2) .* v.^i / sqrt(i!),  i = 0..NumTerms-1
%   (a truncated Taylor feature map of the Gaussian kernel), and the
%   NumTerms^2 Pearson correlations between the features of x and of y are
%   summarised by a p-norm:
%       "max" - largest |r| (IDS-max, the reference default)
%       1     - mean |r|    (IDS-1)
%       2     - root mean r^2 (IDS-2)
%   Cost is NumTerms^2 matrix products of size dx-by-n times n-by-dy, so it
%   scales like CORR by a constant factor.
%
%   Name-value arguments:
%     NumTerms        - Taylor terms per variable. Default 6.
%     PNorm           - "max", 1 or 2. Default "max".
%     Bandwidth       - constant B in exp(-B*v^2). Default 0.5, i.e. a unit
%                       Gaussian kernel.
%     NumPermutations - permutations for empirical p-values. Default 0 (no
%                       p-values). Each column is shuffled independently.
%     ScaleRange      - [lo hi] to min-max scale every column into before the
%                       feature map. Default [] (no scaling, as in the
%                       reference library). Use [0 1] for expression data.
%
%   The kernel bandwidth is fixed and centred at zero, so IDS depends on the
%   origin and range of the input. For scRNA-seq pass library-size
%   normalised, log1p expression with ScaleRange=[0 1] (see NET.IDSNET for
%   the measurements); unscaled log values score near chance on sparse data.
%
%   p-values are (1 + #{null >= observed}) / (1 + NumPermutations), which
%   never returns 0; the reference library uses #{null > observed} / N.
%
%   Ref: Radhakrishnan, Jain, Uhler & Lander, PNAS 122:e2509860122 (2025).
%   Port of https://github.com/aradha/interdependence_scores (numpy backend).
%
%   See also CORR, PKG.E_XICOR, NET.IDSNET.

arguments
    X {mustBeNumeric, mustBeReal}
    Y {mustBeNumeric, mustBeReal} = []
    options.NumTerms (1,1) {mustBeInteger, mustBePositive} = 6
    options.PNorm = "max"
    options.Bandwidth (1,1) {mustBeNumeric, mustBePositive} = 0.5
    options.NumPermutations (1,1) {mustBeInteger, mustBeNonnegative} = 0
    options.ScaleRange {mustBeNumeric} = []
end

pnorm = i_checkpnorm(options.PNorm);
X = i_prepare(X, options.ScaleRange);
hasY = ~isempty(Y);
if hasY
    Y = i_prepare(Y, options.ScaleRange);
    if size(Y, 1) ~= size(X, 1)
        error("X and Y must have the same number of rows (samples).");
    end
end

k = options.NumTerms;
B = options.Bandwidth;
if hasY
    C = i_ids(X, Y, k, B, pnorm);
else
    C = i_ids(X, [], k, B, pnorm);
end

p = [];
if nargout > 1 && options.NumPermutations > 0
    nperm = options.NumPermutations;
    exceed = zeros(size(C));
    for t = 1:nperm
        Xp = i_shufflecolumns(X);
        if hasY
            null = i_ids(Xp, i_shufflecolumns(Y), k, B, pnorm);
        else
            null = i_ids(Xp, [], k, B, pnorm);
        end
        exceed = exceed + (null >= C);
    end
    p = (1 + exceed)./(1 + nperm);
    if ~hasY
        % A variable is trivially dependent on itself, but shuffling it also
        % shuffles its partner, so the permutation null says nothing here.
        p(1:size(p, 1) + 1:end) = 1/(1 + nperm);
    end
elseif nargout > 1
    warning("Set NumPermutations > 0 to compute p-values; returning [].");
end
end


function C = i_ids(X, Y, k, B, pnorm)
% Standardised features: each column centred and scaled to unit norm, so
% Za{a}'*Zb{b} holds the Pearson correlations between term a and term b.
Zx = i_features(X, k, B);
if isempty(Y)
    Zy = Zx;
else
    Zy = i_features(Y, k, B);
end
symmetric = isempty(Y);

C = zeros(size(X, 2), size(Zy{1}, 2), "like", Zx{1});
for a = 1:k
    if symmetric
        bset = a:k;
    else
        bset = 1:k;
    end
    for b = bset
        R = abs(Zx{a}'*Zy{b});
        if symmetric && b > a
            % Term pair (b, a) is the transpose of (a, b).
            R = i_accumulate2(R, R', pnorm);
        end
        switch pnorm
            case "max"
                C = max(C, R);
            case "1"
                C = C + R;
            case "2"
                C = C + R.^2;
            otherwise
                % pnorm is validated by i_checkpnorm; no other case occurs.
        end
    end
end

switch pnorm
    case "1"
        C = C/k^2;
    case "2"
        C = sqrt(C/k^2);
    otherwise
        % "max" needs no normalisation.
end
end


function R = i_accumulate2(R1, R2, pnorm)
% Merge the (a,b) and (b,a) blocks so the caller counts both.
switch pnorm
    case "max"
        R = max(R1, R2);
    case "1"
        R = R1 + R2;
    case "2"
        R = sqrt(R1.^2 + R2.^2);
    otherwise
        R = R1;
end
end


function Z = i_features(V, k, B)
envelope = exp(-B*V.^2);
Z = cell(1, k);
for i = 0:k-1
    F = envelope.*V.^i/sqrt(factorial(i));
    F = F - mean(F, 1);
    s = sqrt(sum(F.^2, 1));
    F = F./s;
    % Constant features (e.g. an all-zero gene) have no defined
    % correlation; the reference maps the resulting NaNs to 0.
    F(:, s <= eps(class(F))) = 0;
    Z{i+1} = F;
end
end


function V = i_prepare(V, scaleRange)
if issparse(V)
    V = full(V);
end
if ~isfloat(V)
    V = double(V);
end
if ~isempty(scaleRange)
    lo = min(V, [], 1);
    span = max(V, [], 1) - lo;
    span(span == 0) = 1;
    V = scaleRange(1) + (scaleRange(2) - scaleRange(1))*(V - lo)./span;
end
end


function V = i_shufflecolumns(V)
[n, d] = size(V);
[~, idx] = sort(rand(n, d), 1);
V = V(idx + (0:d-1)*n);
end


function pnorm = i_checkpnorm(pnorm)
if (isstring(pnorm) || ischar(pnorm)) && strcmpi(pnorm, "max")
    pnorm = "max";
elseif isnumeric(pnorm) && isscalar(pnorm) && ismember(pnorm, [1 2])
    pnorm = string(pnorm);
else
    error("PNorm must be ""max"", 1 or 2.");
end
end
