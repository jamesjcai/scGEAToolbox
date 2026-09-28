function [A, info] = i_cpals(XM, K, opts)
%I_CPALS  Denoised network by R scTenifoldNet's CP-ALS tensor decomposition.
%   A = TEN.I_CPALS(XM, K) decomposes the genes x genes x networks tensor XM
%   into K rank-one components and returns the mean of the estimated
%   slices, divided by its largest absolute value and rounded, as
%   tensorDecomposition() in the R package scTenifoldNet.
%
%   The decomposition is a port of the package's cpDecomposition(): the
%   factors are initialized with R's rnorm after set.seed(Seed) (through
%   TEN.RRANDOM), each column is normalized by its L1 norm, and iteration
%   stops when the change of the residual norm, relative to the norm of the
%   tensor, falls below Tol, or after MaxIter - 1 iterations. Unlike
%   TEN.I_TD1 it needs no Tensor Toolbox.
%
%   [A, INFO] = TEN.I_CPALS(...) also returns a struct with Lambdas, U (the
%   factor matrices), Converged, Residuals and NormPercent.
%
%   NAME-VALUE ARGUMENTS:
%     NumDecimals - digits A is rounded to (default 1; scTenifoldKnk uses 3)
%     MaxIter     - maximum iterations (default 1000)
%     Tol         - relative tolerance (default 1e-5)
%     Seed        - seed of the initial factors, as set.seed (default 1)
%
% see also: TEN.RRANDOM, TEN.I_TD1, TEN.SCTENIFOLDNET

arguments
    XM (:, :, :) double
    K (1, 1) double {mustBeInteger, mustBePositive} = 3
    opts.NumDecimals (1, 1) double {mustBeInteger} = 1
    opts.MaxIter (1, 1) double {mustBePositive} = 1000
    opts.Tol (1, 1) double {mustBePositive} = 1e-5
    opts.Seed (1, 1) double {mustBeInteger} = 1
end

[nI, nJ, nSlices] = size(XM);
tensorNormSq = sum(XM.*XM, "all");
tensorNorm = sqrt(tensorNormSq);

% R fills each factor column by column from one rnorm() call per mode
stream = ten.RRandom(opts.Seed);
U1 = reshape(stream.rnorm(nI*K), nI, K);
U2 = reshape(stream.rnorm(nJ*K), nJ, K);
U3 = reshape(stream.rnorm(nSlices*K), nSlices, K);

iter = 1;
converged = false;
residuals = zeros(opts.MaxIter, 1);
lambdas = zeros(1, K);
prevResid = Inf;
while iter < opts.MaxIter && ~converged
    % Mode 1, slicewise MTTKRP. inv() rather than "/" follows R's solve(V).
    V = (U2.'*U2).*(U3.'*U3);
    mttkrp = zeros(nI, K);
    for k = 1:nSlices
        mttkrp = mttkrp + (XM(:, :, k)*U2).*U3(k, :);
    end
    U1 = i_normalize(mttkrp*inv(V)); %#ok<MINV>

    V = (U1.'*U1).*(U3.'*U3);
    mttkrp = zeros(nJ, K);
    for k = 1:nSlices
        mttkrp = mttkrp + (XM(:, :, k).'*U1).*U3(k, :);
    end
    U2 = i_normalize(mttkrp*inv(V)); %#ok<MINV>

    V = (U1.'*U1).*(U2.'*U2);
    mttkrp = zeros(nSlices, K);
    for k = 1:nSlices
        mttkrp(k, :) = sum(U1.*(XM(:, :, k)*U2), 1);
    end
    [U3, lambdas] = i_normalize(mttkrp*inv(V)); %#ok<MINV>

    % ||X - est||^2 = ||X||^2 - 2<X, est> + ||est||^2, without
    % reconstructing the tensor
    innerProduct = sum(lambdas.*sum(U3.*mttkrp, 1));
    gram = (U1.'*U1).*(U2.'*U2).*(U3.'*U3);
    estNormSq = sum((lambdas.'*lambdas).*gram, "all");
    currResid = sqrt(max(tensorNormSq - 2*innerProduct + estNormSq, 0));

    residuals(iter) = currResid;
    if iter > 1 && abs(currResid - prevResid)/tensorNorm < opts.Tol
        converged = true;
    else
        prevResid = currResid;
        iter = iter + 1;
    end
end

% Sum of the estimated slices, one at a time as in R
A = zeros(nI, nJ);
for k = 1:nSlices
    A = A + (U1.*(lambdas.*U3(k, :)))*U2.';
end
A = A/nSlices;
A = A/max(abs(A), [], "all");
A = round(A, opts.NumDecimals);

if nargout > 1
    residuals = residuals(residuals ~= 0);
    if isempty(residuals)
        normPercent = NaN;
    else
        normPercent = (1 - residuals(end)/tensorNorm)*100;
    end
    info = struct("Lambdas", lambdas, "U", {{U1, U2, U3}}, ...
        "Converged", converged, "Residuals", residuals, ...
        "NormPercent", normPercent);
end
end

function [U, lambdas] = i_normalize(U)
lambdas = sum(abs(U), 1);
U = U./lambdas;
end
