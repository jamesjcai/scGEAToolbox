function [Cp,dp,x,bestx,bestf] = prepareQuboForTabu(obj, C, x)
%PREPAREQUBOFORTABU Prepares a QUBO for a call to Tabu Search

%   [CP, DP, X, BESTX, BESTF] = PREPAREQUBOFORTABU(OBJ, C, X) prepares
%   the Quadratic Unconstrained Binary Problem defined by quadratic term C
%   and linear term d to be called by Tabu search 

%   Copyright 2022-2023 The MathWorks, Inc.


% Step 0: Transform problem into the form tabu search requires
% TODO: Keep the matrix form for now. However, should we just be reusing
% QUBOProblem as the container to avoid memory copies (?)
[Cp, dp] = quantum.internal.optimization.tabusearch.MultiCallTabuSearch.moveDiagonal2Linear(C);

% Density of the quadratic matrix
density = nnz(C)/numel(C);
if obj.StoreQSparse && density < 0.25
    Cp = quantum.internal.optimization.tabusearch.SparseFixedPositionSymmetric(Cp);
else
    Cp = quantum.internal.optimization.tabusearch.DenseFixedPositionSymmetric(full(Cp));
end

% Step 1: Initial x is passed

bestx = x;
bestf = x'*C*x;

% Step 2: Transform the problem instance according to formulas
% (3), (4) in Paulbeckis (2006) applied with respect to x.
[Cp, dp] = quantum.internal.optimization.tabusearch.MultiCallTabuSearch.transformProblem(Cp, dp, x);

end

