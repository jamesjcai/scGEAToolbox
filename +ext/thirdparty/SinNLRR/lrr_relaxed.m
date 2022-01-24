function [ Z,it_num] = lrr_relaxed( X,lambda,mu)
% min || X - XZ ||_F^2 + lambda || Z ||_*
% via ADMM
%

if (nargin<3)
    mu = 1;
end

max_iterations = 200;

%func_vals = zeros(max_iterations, 1);

n = size(X, 2);

Z = zeros(n);
J = zeros(n);
Y = zeros(n);

var_X=X'*X;

A=inv(var_X + mu*speye((n)));

tol_1 = 1*10^-4;

for k = 1 : max_iterations
    
    Z_prev = Z;
    J_prev = J;
        
    J = A*(var_X + Y + mu*Z);
    J(J<0)=0;
    
    V = J - (1/mu) * Y;
    
    [Z, ~] = solve_nn(V, lambda / mu );
    Z(Z<0)=0;
    Y = Y + mu*(Z - J);
    
    % Check convergence
    if (max(max(abs(J - Z))) < tol_1&&max(max(abs(Z - Z_prev))) < tol_1&&max(max(abs(J - J_prev))) < tol_1)
        break;
    end
    
end
it_num=k;
end