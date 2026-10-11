function [lambda, U, info] = i_cpalstt(XM, R, opts)
%I_CPALSTT  Tensor Toolbox's cp_als, on a stack of networks, without the toolbox.
%   [LAMBDA, U] = TEN.I_CPALSTT(XM, R) decomposes the n1 x n2 x K array XM
%   into R rank-one components exactly as Tensor Toolbox 3.8's
%
%       cp_als(tensor(XM), R, 'init', 'nvecs', 'maxiters', 100)
%
%   does, and returns its weights LAMBDA (R x 1) and factor matrices
%   U = {U1, U2, U3}, so that slice k of the estimate is
%   U1*diag(LAMBDA.*U3(k, :)')*U2'. This is the decomposition TEN.DO_TD_CP
%   and so scTenifoldNet, scTenifoldKnk and net.pcrnet_denoised use.
%
%   [LAMBDA, U, INFO] = TEN.I_CPALSTT(...) also returns a struct with
%   Iterations, Fit and Converged.
%
%   NAME-VALUE ARGUMENTS:
%     MaxIter - maximum iterations (default 100, as TEN.DO_TD_CP set it)
%     Tol     - stop when the fit changes by less than this (default 1e-4)
%     Verbose - print the fit at every iteration, as cp_als's printitn = 1
%               (default true)
%
%   WHY A PORT. cp_als works on the dense tensor: every MTTKRP permutes and
%   reshapes all n1*n2*K entries, and nvecs forms the n x nK unfolding. The
%   networks scTenifoldNet stacks are about 95% zeros after its q = 0.95
%   filter, so here the same products are taken on one sparse unfolding of
%   the stack, nvecs's Gram matrices are built without the dense unfolding,
%   and only the mean slice need ever be formed (TEN.I_TD1). At 3000 genes
%   and 10 networks: 5.6 s down to about 1.2 s.
%   The arithmetic is cp_als's, step for step: the nvecs start on modes 2
%   and 3 (eigs, largest magnitude, each vector's largest entry made
%   positive), the 2-norm on the first iteration and the max-norm (at least
%   1) after it, and the same fit and stopping rule. cp_als's closing ARRANGE
%   and FIXSIGNS only reorder components and flip signs in pairs, which does
%   not change the estimate, so they are left out. Results agree with
%   cp_als to rounding error; tests/cpalsttTest.m compares the two.
%
% See also TEN.DO_TD_CP, TEN.I_TD1, TEN.I_CPALS

arguments
    XM (:, :, :) {mustBeNumeric}
    R (1, 1) double {mustBeInteger, mustBePositive}
    opts.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    opts.Tol (1, 1) double {mustBePositive} = 1e-4
    opts.Verbose (1, 1) logical = true
end

[n1, n2, K] = size(XM);
XM = double(XM);

% The stack unfolded once, H = [X_1, ..., X_K] (n1 x n2*K), sparse when
% the networks are mostly zeros, as filtered PCNet stacks are, and dense
% otherwise, where sparse storage would only cost. Each MTTKRP below is then
% one product with H, not K products with the slices, and mode 3
% reuses mode 2's product (2026-10-05: 2.0 s -> about 1.2 s at 3000 genes).
useSparse = nnz(XM) < 0.25*numel(XM);
if useSparse
    H = sparse(reshape(XM, n1, n2*K));
else
    H = reshape(XM, n1, n2*K);
end

% Mode-3 Gram matrix, the slices' inner products, in one dense product;
% its trace is the tensor's squared norm.
V = reshape(XM, n1*n2, K);
G3 = V.'*V;
G3 = (G3 + G3.')/2;
clear V
normXsq = trace(G3);
normX = sqrt(normXsq);

% nvecs on modes 2 and 3 (dimorder(2:end)); mode 1 is solved first.
% Mode-2 Gram matrix: sum_k X_k'*X_k, dense even when the slices are not;
% eigs on the matrix (as nvecs calls it) beat applying the stack inside
% eigs, 0.6 s against 0.9 s at 3000 genes.
G2 = zeros(n2, n2);
for k = 1:K
    G2 = G2 + full(H(:, (k - 1)*n2 + (1:n2)).'*XM(:, :, k));
end
clear XM
[U2, ~] = eigs((G2 + G2.')/2, R, "largestabs", "Display", false);
clear G2
U2 = i_flipsign(U2);
[U3, ~] = eigs(G3, R, "largestabs", "Display", false);
U3 = i_flipsign(U3);
U1 = zeros(n1, R);

UtU = zeros(R, R, 3);
UtU(:, :, 2) = U2.'*U2;
UtU(:, :, 3) = U3.'*U3;

if opts.Verbose
    storage = "dense";
    if useSparse
        storage = "sparse";
    end
    fprintf('CP-ALS: tensor %s, %s slices, R = %d, tol = %g\n', ...
        mat2str([n1, n2, K]), storage, R, opts.Tol);
end
fit = 0;
converged = false;
for iter = 1:opts.MaxIter
    fitold = fit;

    % Mode 1: sum_k X_k*U2 .* U3(k, :) = H * (U3 Khatri-Rao U2)
    KR = reshape(reshape(U2, n2, 1, R).*reshape(U3, 1, K, R), n2*K, R);
    M = full(H*KR);
    [U1, ~] = i_update(M, UtU(:, :, 2).*UtU(:, :, 3), iter);
    UtU(:, :, 1) = U1.'*U1;

    % Mode 2: sum_k X_k'*U1 .* U3(k, :). Q(:, k, :) holds X_k'*U1.
    Q = reshape(full((U1.'*H).'), n2, K, R);
    M = reshape(sum(Q.*reshape(U3, 1, K, R), 2), n2, R);
    [U2, ~] = i_update(M, UtU(:, :, 1).*UtU(:, :, 3), iter);
    UtU(:, :, 2) = U2.'*U2;

    % Mode 3: row k is sum(U1 .* (X_k*U2)) = sum(U2 .* (X_k'*U1)), from Q
    M3 = reshape(sum(Q.*reshape(U2, n2, 1, R), 1), K, R);
    [U3, lambda] = i_update(M3, UtU(:, :, 1).*UtU(:, :, 2), iter);
    UtU(:, :, 3) = U3.'*U3;

    iprod = sum(U3.*M3, 1)*lambda;
    normPsqr = sum(UtU(:, :, 1).*UtU(:, :, 2).*UtU(:, :, 3).*(lambda*lambda.'), "all");
    if normX == 0
        fit = normPsqr - 2*iprod;
    else
        fit = 1 - sqrt(normXsq + normPsqr - 2*iprod)/normX;
    end
    fitchange = abs(fitold - fit);
    converged = iter > 1 && fitchange < opts.Tol;
    if opts.Verbose
        fprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange);
    end
    if converged
        break
    end
end

U = {U1, U2, U3};
info = struct("Iterations", iter, "Fit", fit, "Converged", converged);
end

function u = i_flipsign(u)
% nvecs's flipsign: the largest-magnitude entry of each vector positive
[~, loc] = max(abs(u), [], 1);
idx = sub2ind(size(u), loc, 1:size(u, 2));
u(:, u(idx) < 0) = -u(:, u(idx) < 0);
end

function [Unew, lambda] = i_update(M, Y, iter)
Unew = M/Y;
if iter == 1
    lambda = sqrt(sum(Unew.^2, 1)).';   % 2-norm
else
    lambda = max(max(abs(Unew), [], 1), 1).';   % max-norm
end
Unew = Unew./lambda.';
end
