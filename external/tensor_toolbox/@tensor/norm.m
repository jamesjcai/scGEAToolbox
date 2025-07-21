function n = norm(T)
%NORM Frobenius norm of a tensor.
%
%   NORM(X) returns the Frobenius norm of a tensor.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


v = reshape(T.data, numel(T.data), 1);
n = norm(v);
