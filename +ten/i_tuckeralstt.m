function [core, U, info] = i_tuckeralstt(XM, R, opts)
%I_TUCKERALSTT  Tensor Toolbox's tucker_als, on a stack of networks, without the toolbox.
%   [CORE, U] = TEN.I_TUCKERALSTT(XM, R) computes the rank-(R, R, R) Tucker
%   decomposition of the n1 x n2 x K array XM exactly as Tensor Toolbox
%   3.8's
%
%       tucker_als(tensor(XM), R, 'printitn', 0)
%
%   does (higher-order orthogonal iteration from a random start), and
%   returns the R x R x R core and the factor matrices U = {U1, U2, U3}, so
%   that slice k of the estimate is U1*(sum_c CORE(:, :, c)*U3(k, c))*U2'.
%   This is the 'Tucker' option of scTenifoldNet, scTenifoldKnk and
%   net.pcrnet_denoised, through TEN.DO_TD_TUCKER.
%
%   [CORE, U, INFO] = TEN.I_TUCKERALSTT(...) also returns a struct with
%   Iterations and Fit.
%
%   NAME-VALUE ARGUMENTS:
%     MaxIter - maximum iterations (default 50, tucker_als's)
%     Tol     - stop when the fit changes by less than this (default 1e-4)
%
%   As in tucker_als, the start is rand(n2, R) and then rand(K, R) from
%   the global random stream, so the same seed gives the same result. Each
%   factor is the leading R left singular vectors of the matching unfolding
%   of the partly projected tensor - nvecs, with each vector's
%   largest-magnitude entry made positive - computed on one sparse
%   unfolding of the stack instead of on the dense tensor.
%   tests/tensorPortsTest.m compares the two.
%
% See also TEN.DO_TD_TUCKER, TEN.I_CPALSTT, TEN.I_TD1

arguments
    XM (:, :, :) {mustBeNumeric}
    R (1, 1) double {mustBeInteger, mustBePositive}
    opts.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 50
    opts.Tol (1, 1) double {mustBePositive} = 1e-4
end

[n1, n2, K] = size(XM);
XM = double(XM);
% The stack laid out twice, side by side, H = [X_1, ..., X_K] (n1 x n2*K),
% and stacked, B = [X_1; ...; X_K] (n1*K x n2): sparse when the networks
% are mostly zeros and dense otherwise. Mode 1 is then one product with B
% and mode 2 one with H, rather than K products with the slices each, and
% mode 3 reuses mode 2's product (2026-10-05).
useSparse = nnz(XM) < 0.25*numel(XM);
if useSparse
    H = sparse(reshape(XM, n1, n2*K));
    [i, jk, v] = find(H);
    j = mod(jk - 1, n2) + 1;
    B = sparse(i + n1*(jk - j)/n2, j, v, n1*K, n2);
    clear i jk v j
else
    H = reshape(XM, n1, n2*K);
    B = reshape(permute(XM, [1 3 2]), n1*K, n2);
end
clear XM
normXsq = full(sum(H.^2, "all"));

% tucker_als's 'random' start, modes 2 and 3 in that order
U2 = rand(n2, R);
U3 = rand(K, R);

fit = 0;
for iter = 1:opts.MaxIter
    fitold = fit;

    % Mode 1: X x2 U2' x3 U3', unfolded n1 x R^2 (mode-2 index fastest).
    % P(i, k, a) is (X_k*U2)(i, a); column a + R*(b-1) of the unfolding is
    % sum_k P(:, k, a)*U3(k, b).
    P = reshape(full(B*U2), n1, K, R);
    T1 = reshape(reshape(permute(P, [1 3 2]), n1*R, K)*U3, n1, R*R);
    U1 = i_nvecs(T1, R);

    % Mode 2: X x1 U1' x3 U3', unfolded n2 x R^2 (mode-1 index fastest).
    % Q(j, k, a) is (X_k'*U1)(j, a).
    Q = reshape(full((U1.'*H).'), n2, K, R);
    T2 = reshape(reshape(permute(Q, [1 3 2]), n2*R, K)*U3, n2, R*R);
    U2 = i_nvecs(T2, R);

    % Mode 3: X x1 U1' x2 U2', unfolded K x R^2. W(a, b, k) is
    % U1(:, a)'*X_k*U2(:, b) = sum_j Q(j, k, a)*U2(j, b).
    W = permute(reshape(U2.'*reshape(Q, n2, K*R), R, K, R), [3 1 2]);
    T3 = reshape(W, R*R, K).';
    U3 = i_nvecs(T3, R);

    core = reshape(reshape(W, R*R, K)*U3, R, R, R);
    normresidual = sqrt(max(normXsq - sum(core.^2, "all"), 0));
    fit = 1 - normresidual/sqrt(normXsq);
    if iter > 1 && abs(fitold - fit) < opts.Tol
        break
    end
end

U = {U1, U2, U3};
info = struct("Iterations", iter, "Fit", fit);
end

function u = i_nvecs(Tn, r)
% nvecs: leading eigenvectors of Tn*Tn', i.e. left singular vectors of the
% thin unfolding, with each vector's largest-magnitude entry positive
[u, ~, ~] = svd(Tn, "econ");
u = u(:, 1:r);
[~, loc] = max(abs(u), [], 1);
idx = sub2ind(size(u), loc, 1:r);
u(:, u(idx) < 0) = -u(:, u(idx) < 0);
end
