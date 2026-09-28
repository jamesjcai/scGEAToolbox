function XM = i_ncref(X, opts)
%I_NCREF  Subsampled PC networks, as makeNetworks() in R scTenifoldNet.
%   XM = TEN.I_NCREF(X) builds NumNets principal-component regression
%   networks from the genes-by-cells matrix X, each on NumCells cells drawn
%   with replacement, and returns them as a genes x genes x NumNets array.
%   With the same seed it draws the same cells as
%
%       set.seed(Seed); makeNetworks(X, nNet, nCells, nComp, q = Q)
%
%   and gives the same networks up to floating-point rounding. It differs
%   from TEN.I_NC, the default builder, in four ways, all following R:
%     - cells are drawn by TEN.RRANDOM, i.e. R's sample(), always with
%       replacement (no jackknife, no fallback on small populations)
%     - genes not expressed in the drawn cells are left out of the
%       regression and get no edges
%     - edges are filtered by R's type-7 quantile of |A| over the block of
%       non-constant genes, diagonal included, keeping |A| >= threshold
%     - scaling is by max(|A|) (TEN.E_FILTADJC uses the top-10 mean)
%
%   NAME-VALUE ARGUMENTS:
%     NumNets  - number of networks (default 10)
%     NumCells - cells per network (default 500)
%     NumComp  - principal components per regression (default 3)
%     Q        - edge-filtering quantile; 0 or 1 disable (default 0.95)
%     Seed     - seed, as set.seed (default 1)
%
% see also: TEN.I_NC, TEN.RRANDOM, NET.PCRNET, TEN.SCTENIFOLDNET

arguments
    X double
    opts.NumNets (1, 1) double {mustBeInteger, mustBePositive} = 10
    opts.NumCells (1, 1) double {mustBeInteger, mustBePositive} = 500
    opts.NumComp (1, 1) double {mustBeInteger} = 3
    opts.Q (1, 1) double {mustBeBetween(opts.Q, 0, 1)} = 0.95
    opts.Seed (1, 1) double {mustBeInteger} = 1
end

[nGenes, nCells] = size(X);
if opts.NumComp < 2 || opts.NumComp >= nGenes
    error("ten:i_ncref:badNumComp", ...
        "NumComp should be >= 2 and < the number of genes (%d). Received %d.", ...
        nGenes, opts.NumComp);
end

stream = ten.RRandom(opts.Seed);
XM = zeros(nGenes, nGenes, opts.NumNets);
for k = 1:opts.NumNets
    fprintf('Building network...%d of %d\n', k, opts.NumNets);
    Z = X(:, stream.sample(nCells, opts.NumCells));
    expressed = sum(Z, 2) > 0;
    Z = Z(expressed, :);

    % pcNet: regression, then scale, then filter. net.pcrnet returns zero
    % rows and columns for constant genes; pcNet's quantile is over the
    % non-constant ones only.
    A = net.pcrnet(Z, opts.NumComp);
    used = any(Z ~= Z(:, 1), 2);
    B = A(used, used);
    mx = max(abs(B), [], "all");
    if isfinite(mx) && mx > 0
        B = B/mx;
    end
    if opts.Q > 0 && opts.Q < 1
        absB = abs(B);
        B(absB < i_quantile7(absB(:), opts.Q)) = 0;
    end
    A(used, used) = B;
    XM(expressed, expressed, k) = A;
end
end

function qs = i_quantile7(x, p)
% R's quantile(x, p), type 7, with R's interpolation formula
n = numel(x);
x = sort(x);
index = 1 + max(n - 1, 0)*p;
lo = floor(index);
hi = ceil(index);
qs = x(lo);
if index > lo && x(hi) ~= qs
    h = index - lo;
    qs = (1 - h)*qs + h*x(hi);
end
end
