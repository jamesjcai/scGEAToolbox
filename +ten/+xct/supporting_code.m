% https://github.com/cailab-tamu/scTenifoldXct/blob/main/scTenifoldXct/stiefel.py


function S = supporting_code(M, Z)
% Return the skew-symmetric part of M'*Z.
S = 0.5 * (M'*Z - Z'*M);
end

% function P = proj_stiefel(M, Z)
% % Project Z onto the Stiefel manifold defined by M.
% MskewMTZ = M * skew(M, Z);
% IMMTZ = (eye(size(M)) - M * M') * Z;
% P = MskewMTZ + IMMTZ;
% end

% function Q = rand_stiefel(n, p)
% % Generate a random Stiefel point using QR decomposition of a random Gaussian matrix.
% X = randn(n, p);
% [Q, ~] = qr(X, 0);  % Ensure economy-size QR decomposition
% end

% function P = retr_stiefel(Z)
% % Retract Z onto Stiefel manifold using truncated SVD.
% [U, ~, V] = svd(Z, 'econ');
% P = U * V';
% end


