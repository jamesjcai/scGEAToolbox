function Xhat0 = do_td_tucker(XM0)
% Rank-(5, 5, 5) Tucker decomposition of a stack of networks (n x n x k);
% returns the estimated stack. ten.i_tuckeralstt is a port of Tensor
% Toolbox's tucker_als(tensor(XM0), 5, 'printitn', 0): same random start
% from the global stream, same result to rounding error.

[core, U] = ten.i_tuckeralstt(XM0, 5);
K = size(XM0, 3);
Xhat0 = zeros(size(U{1}, 1), size(U{2}, 1), K);
for k = 1:K
    Xhat0(:, :, k) = U{1}*sum(core.*reshape(U{3}(k, :), 1, 1, []), 3)*U{2}.';
end
end
