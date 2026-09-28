function [A] = pcrnet(X, ncom, fastersvd, dozscore, UseParallel, guiwaitbar, UseGPU)
% Construct GRN using principal component regression (PCR)
%
% A = net.pcrnet(X)
% A = net.pcrnet(X, ncom)
% A = net.pcrnet(X, ncom, fastersvd, dozscore, UseParallel, guiwaitbar, UseGPU)
%
% X           - genes x cells expression matrix (LogNormalized recommended)
% ncom        - number of principal components (default: 3)
% fastersvd   - ignored; kept for compatibility (default: false)
% dozscore    - z-score genes before regression (default: true)
% UseParallel - ignored; kept for compatibility (default: false)
% guiwaitbar  - show GUI progress bar (default: false)
% UseGPU      - use CUDA GPU via gpuArray (default: false)
%
% Every gene k is regressed on the top ncom principal components of all the
% other genes; row k of A holds the coefficients. Genes that are constant
% across cells (after z-scoring) get zero rows and columns.
%
% Rather than one SVD per gene, all regressions come from one
% eigendecomposition. With G = X*X' = Q*diag(d)*Q', removing gene y gives
% G - y*y', i.e. diag(d) - z*z' in the eigenbasis (z = Q'*y). Its eigenvalues
% are the roots of the secular equation 1 - sum(z.^2 ./ (d - mu)) = 0, one in
% each interval (d(i+1), d(i)), with eigenvectors proportional to z ./ (d - mu).
% The result is exact, not an approximation of the per-gene SVD.
% Port of pcNet() in scTenifoldNet 1.4.1 (Daniel Osorio, 2026).
%
% ref: https://github.com/cailab-tamu/scTenifoldNet/blob/master/R/pcNet.R
%      https://rdrr.io/cran/dna/man/PCnet.html

arguments
    X double
    ncom(1, 1) {mustBeNumeric} = 3
    fastersvd(1, 1) logical = false %#ok<INUSA>
    dozscore(1, 1) logical = true
    UseParallel(1, 1) logical = false %#ok<INUSA>
    guiwaitbar(1, 1) logical = false
    UseGPU(1, 1) logical = false
end

if UseGPU && gpuDeviceCount < 1
    warning('net:pcrnet:NoGPU', 'No CUDA GPU found; falling back to CPU.');
    UseGPU = false;
end

% cells x genes, dense
X = full(X.');
nGenes = size(X, 2);
nCells = size(X, 1);

% Genes without information are set aside: constant genes when z-scoring
% (they z-score to zero), all-zero genes otherwise. Tested exactly, before
% standardizing, so rounding in the mean cannot turn a constant into noise.
if dozscore
    used = any(X ~= X(1, :), 1);
else
    used = any(X ~= 0, 1);
end
X = X(:, used);
nUsed = size(X, 2);

if ncom >= nUsed
    error('net:pcrnet:TooFewGenes', ...
        'ncom (%d) must be smaller than the number of informative genes (%d).', ncom, nUsed);
end

if dozscore
    X = X - mean(X, 1);
    X = X ./ sqrt(sum(X.*X, 1)/max(1, nCells - 1));
end
if UseGPU
    X = gpuArray(X);
end

if guiwaitbar
    fw = gui.myWaitbar([]);
    progress = @(f) gui.myWaitbar([], fw, [], [], '', f);
else
    progress = @(f) [];
end

% Z: r x genes, gene j expressed in the eigenbasis of the Gram matrix;
% d: eigenvalues in decreasing order. The smaller Gram matrix is used.
if nCells <= nUsed
    G = X*X.';
    [Q, d] = eig((G + G.')/2, 'vector');
    [d, order] = sort(max(d, 0), 'descend');
    Z = Q(:, order).'*X;
else
    G = X.'*X;
    [V, d] = eig((G + G.')/2, 'vector');
    [d, order] = sort(max(d, 0), 'descend');
    Z = sqrt(d).*V(:, order).';
end
clear X G Q V

W = secularWeights(Z, d, ncom, progress);
B = W.'*Z;
B(1:(nUsed + 1):end) = 0;

A = zeros(nGenes, nGenes);
A(used, used) = gather(B);
if guiwaitbar
    gui.myWaitbar([], fw);
end
end

function W = secularWeights(Z, d, ncom, progress)
% Solves the leave-one-gene-out eigenproblems for all genes at once.
%
% Column k of Z is z = Q'*y for gene k. Returns W (r x genes) with column k
% equal to sum_i u_i*(u_i'*z)/mu_i over the top ncom eigenpairs (mu_i, u_i)
% of diag(d) - z*z'; the coefficients of gene k are then Z'*W(:, k).
%
% Each root mu_i is the only root in (d(i+1), d(i)) of the decreasing
% f(mu) = 1 - sum(z.^2 ./ (d - mu)). Following LAPACK dlaed4, mu is written
% as origin + tau with origin the nearer interval end, so d - mu is computed
% without cancellation; each step models the sums above and below the
% interval by one-pole functions and solves the resulting quadratic,
% falling back to bisection inside a maintained bracket. With v = z./(d - mu),
% f(mu) = 0 gives u_i'*z = 1/norm(v), so each term of W is v/(mu*norm(v)^2).

maxIterations = 300;
tol = 4*eps;
nEigInput = size(Z, 1);
nGenes = size(Z, 2);

% With very few cells, every interval needs a lower end: the remaining
% eigenvalues of the Gram matrix are zero.
nEig = max(nEigInput, ncom + 1);
if nEig > nEigInput
    Z = [Z; zeros(nEig - nEigInput, nGenes, 'like', Z)];
    d = [d; zeros(nEig - nEigInput, 1, 'like', d)];
end

% A zero z_j leaves d_j an eigenvalue whose eigenvector contributes nothing.
% Raising such entries to a negligible size gives the same result through
% the general formula.
tiny = sqrt(max(d(1), realmin))*1e-60;
nearZero = abs(Z) < tiny;
Z(nearZero) = tiny*(1 - 2*(Z(nearZero) < 0));
Z2 = Z.*Z;

W = zeros(nEig, nGenes, 'like', Z);
for i = 1:ncom
    dUpper = d(i);
    dLower = d(i + 1);
    gap = dUpper - dLower;
    % Equal eigenvalues: mu_i equals them and its eigenvector is orthogonal
    % to z, so it contributes nothing.
    if gap > 0
        above = 1:i;
        below = (i + 1):nEig;

        % f(midpoint) >= 0 means the root lies in the upper half
        nearUpper = (1 - sum(Z2./(d - (dLower + gap/2)), 1)) >= 0;
        origin = dLower + gap*nearUpper;

        % Everything is relative to origin (tau = mu - origin), per gene
        poleLower = -gap*nearUpper;
        poleUpper = gap*(~nearUpper);
        bracketLo = -gap/2*nearUpper;
        bracketHi = gap/2*(~nearUpper);
        tau = (bracketLo + bracketHi)/2;
        dMinusOrigin = d - origin;

        active = 1:nGenes;
        for iteration = 1:maxIterations
            tauA = tau(active);
            gaps = dMinusOrigin(:, active) - tauA;
            terms = Z2(:, active)./gaps;
            slopes = terms./gaps;
            sumAbove = sum(terms(above, :), 1);
            slopeAbove = sum(slopes(above, :), 1);
            sumBelow = sum(terms(below, :), 1);
            slopeBelow = sum(slopes(below, :), 1);
            f = 1 - sumAbove - sumBelow;

            % f is decreasing: f > 0 means the root is above tau
            lo = bracketLo(active);
            hi = bracketHi(active);
            lo(f > 0) = tauA(f > 0);
            hi(f < 0) = tauA(f < 0);
            bracketLo(active) = lo;
            bracketHi(active) = hi;

            % One-pole models matching value and slope at tau; solving
            % model = 1 for tau + step gives a*step^2 - b*step + c = 0
            toLower = poleLower(active) - tauA;
            toUpper = poleUpper(active) - tauA;
            a = f + slopeBelow.*toLower + slopeAbove.*toUpper;
            b = a.*(toLower + toUpper) - slopeBelow.*toLower.^2 - slopeAbove.*toUpper.^2;
            c = f.*toLower.*toUpper;

            % Stable quadratic roots; keep the one between the poles
            sqrtDisc = sqrt(max(b.*b - 4*a.*c, 0));
            half = (b + sqrtDisc.*(1 - 2*(b < 0)))/2;
            root1 = half./a;
            root2 = c./half;
            step = root2;
            useRoot1 = root1 > toLower & root1 < toUpper;
            step(useRoot1) = root1(useRoot1);

            tauNew = tauA + step;
            outside = isnan(tauNew) | tauNew <= lo | tauNew >= hi;
            tauNew(outside) = (lo(outside) + hi(outside))/2;

            exact = f == 0;
            tauNew(exact) = tauA(exact);
            converged = exact | abs(tauNew - tauA) <= tol*abs(tauNew) | ...
                (hi - lo) <= tol*max(abs(lo), abs(hi));

            tau(active) = tauNew;
            active = active(~converged);
            if isempty(active)
                break
            end
        end

        mu = origin + tau;
        v = Z./(dMinusOrigin - tau);
        W = W + v./(mu.*sum(v.*v, 1));
    end
    progress(i/(ncom + 1));
end
W = W(1:nEigInput, :);
end
