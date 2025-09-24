function qprob = knapsack2qubo(v, w, M, nvp)
%KNAPSACK2QUBO Convert Knapsack problem to QUBO
%
%   QPROB = KNAPSACK2QUBO(v, w, M) converts the Knapsack problem defined
%   below into a QUBO.
%
%            max v'*x    subject to:   w'*x <= M
%             x                        x(i) binary
%
%   QPROB = KNAPSACK2QUBO(___, Penalty=VALUE) specifies the scalar penalty
%   factor for the linear constraints that are incorporated into the QUBO
%   objective.
% 
%   The Knapsack problem is converted to QUBO following section 5.2 of [1].
%
%   Reference: 
%   [1] Lucas, A., Ising formulations of many NP problems. Frontiers in
%   Physics, 2, 5 (2014). https://arxiv.org/abs/1302.5843

%   Copyright 2025 The MathWorks, Inc.

arguments
    v {quantum.internal.optimization.validators.mustBeDouble, mustBeFinite, mustBeReal, mustBeNonempty} 
    w {quantum.internal.optimization.validators.mustBeDouble, mustBeFinite, mustBePositive, mustBeInteger, quantum.internal.optimization.validators.mustHaveSameNumElements(v, w)}
    M (1, 1) {quantum.internal.optimization.validators.mustBeDouble, mustBeFinite, mustBePositive, mustBeInteger}
    nvp.Penalty (1, 1) {quantum.internal.optimization.validators.mustBeDouble, mustBeFinite, mustBePositive} = 1 
end

v = v(:);
w = w(:);
numVars = numel(v);

% The objective, including the penalized cost is
%
%             N                    N
% F(x, Y) =  sum   -v_i*x_i + P* (sum  w_i*x_i + Y - M )^2    [1]
%           i = 1                i = 1
%
% where Y is an integer slack variable.
% 
% As QUBO requires binary variables, Y is written as
%
%      K
% Y = sum  2^k*y_k
%    k = 0
%
% where y_k is binary and K = floor(log2(M)).

% Create slack coefficients
numSlacks = floor(log2(M)) + 1;
slackCoeffs = 2.^(0:numSlacks-1)';

% The variables in the problem are
% [x_1, x_2, ..., x_N, y_0, y_1, ..., y_K]
%
% Equation [1] can be rewritten as
%
% F(x, y) = Q(x, y) + c(x, y) + d
%
% where
%
%                N                      K
% Q(x, y) = P*( sum  w_i*x_i )^2 + P*( sum  2^k*y_k )^2 + 
%              i = 1                  k = 0
%
%                  N                K
%           2*P*( sum  w_i*x_i )*( sum  2^k*y_k )
%                i = 1            k = 0
%
%            N                       N                        K           
% c(x, y) = sum  -v_i*x_i - 2*P*M*( sum  w_i*x_i ) - 2*P*M*( sum  2^k*y_k )
%          i = 1                   i = 1                    k = 0
%
%
% d = P*M^2

% Quadratic term
v1 = [w; zeros(numSlacks, 1)];
v2 = [zeros(numVars, 1);slackCoeffs];
Q = v1*v2' + v2*v1';
Q(1:numVars, 1:numVars) = w*w';
Q(numVars+1:end, numVars+1:end) = slackCoeffs*slackCoeffs';
Q = Q*nvp.Penalty;

% Linear term
c = [-v - 2*nvp.Penalty*M*w;-2*nvp.Penalty*M*slackCoeffs];

% Constant term
d = nvp.Penalty*M^2;

% Create QUBO
qprob = qubo(Q, c, d);

end
