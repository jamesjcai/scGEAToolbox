function C = mtimes(A, B)
%MTIMES Implement scalar multiplication for a ttensor.
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if ~isa(B, 'ttensor') && isscalar(B)
    C = ttensor(B*A.core, A.u);
elseif ~isa(A, 'ttensor') && isscalar(A)
    C = ttensor(A*B.core, B.u);
else
    error('Use mtimes(full(A),full(B)).');
end
