function Z = minus(X, Y)
%MINUS Binary subtraction (-) for tenmat.
%
%   See also TENMAT, TENMAT/TENMATFUN.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


fun = @minus;

% One argument is a scalar
if ((isscalar(X) || isscalar(Y)))
    if (isscalar(Y)) && isa(X, 'tenmat')
        Z = X;
        Z.data = fun(Z.data, Y);
    else
        Z = Y;
        Z.data = fun(X, Z.data);
    end
    return;
end


% Both arguments are tenmats
Z = tenmat(Y);
if ~(isequal(size(Y), size(Z)))
    error('Tenmat size mismatch.')
end
Z.data = fun(X.data, Z.data);
