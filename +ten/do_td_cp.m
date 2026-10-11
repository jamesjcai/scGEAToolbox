function Xhat0 = do_td_cp(XM0)
% CP decomposition of a stack of networks (n x n x k), one component per
% network, with the leading mode-n singular vectors (nvecs) as the initial
% guess; returns the estimated stack. ten.i_cpalstt is a port of Tensor
% Toolbox's cp_als that agrees with it to rounding error but works slice by
% slice on the mostly-zero networks.

K = size(XM0, 3);
[lambda, U] = ten.i_cpalstt(XM0, K);
Xhat0 = zeros(size(U{1}, 1), size(U{2}, 1), K);
for k = 1:K
    Xhat0(:, :, k) = (U{1}.*(lambda.'.*U{3}(k, :)))*U{2}.';
end
end
