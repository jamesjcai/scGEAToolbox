function x = applyMatToDim(x, M, dim)
% Applies matrix M to dimension dim of x, assuming x is a
% 2-by-2-by-2-by-... array.
%
% Example (apply to a vector):
% x = randn(2,2,2,2,2);
% M = randn(2);
% I = eye(2);
% norm(quantum.internal.gate.krond(I, I, I, I, M) * x(:) - applyMatToDim(x(:), M, 1))
% norm(quantum.internal.gate.krond(I, I, I, M, I) * x(:) - applyMatToDim(x(:), M, 2))
% norm(quantum.internal.gate.krond(I, I, M, I, I) * x(:) - applyMatToDim(x(:), M, 3))
% norm(quantum.internal.gate.krond(I, M, I, I, I) * x(:) - applyMatToDim(x(:), M, 4))
% norm(quantum.internal.gate.krond(M, I, I, I, I) * x(:) - applyMatToDim(x(:), M, 5))
%
% Example (apply to a whole operator matrix):
% applyMatToDim(eye(2^4), [0 1; 1 0], 3)

%   Copyright 2021-2022 The MathWorks, Inc.

szx = size(x);

if dim <= 4
    x = reshape(x, 2^dim, []);
    Mk = kron(M, eye(2^(dim-1)));
    x = Mk*x;
else
    x = reshape(x, 2^(dim-1), 2, []);
    x = pagemtimes(x, 'none', M, 'transpose');
end

x = reshape(x, szx);
