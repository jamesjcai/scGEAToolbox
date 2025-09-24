function [items, totalValue, totalWeight] = quboResult2knapsack(result, v, w, M)
%QUBORESULT2KNAPSACK Convert QUBO Result to Knapsack problem
%
%   [ITEMS, TOTALVALUE, TOTALWEIGHT] = QUBORESULT2KNAPSACK(RESULT, V, W, M)
%   converts the QUBO result to a solution for the Knapsack problem. A
%   numeric index vector of the items placed in the Knapsack, their total
%   value and weight is returned.

%   Copyright 2025 The MathWorks, Inc.

arguments
    result (1, 1) {mustBeA(result, "quboResult")}    
    v {mustBeNumeric, mustBeFinite, mustBeReal, mustBeNonempty} 
    w {mustBeNumeric, mustBeFinite, mustBePositive, mustBeInteger, quantum.internal.optimization.validators.mustHaveSameNumElements(v, w)}
    M (1, 1) {mustBeNumeric, mustBeFinite, mustBePositive, mustBeInteger}
end

% If there is no best solution, just return empty items
if isempty(result.BestX)
    items = [];
    totalValue = [];
    totalWeight = [];
    return
end

% Call outside of arguments block to avoid lint for unused M.
mustBeForProblem(result, v, M);

numVars = numel(v);
items = find(result.BestX(1:numVars));
totalValue = sum(v(items));
totalWeight = sum(w(items));

end

function mustBeForProblem(result, v, M)

numSlackBits = ceil(log2(M));
numItems = numel(v);

if ~isequal(numel(result.BestX), numItems + numSlackBits)
    error(message("quantum:annealing:quboResult2knapsack:mustBeForProblem"));
end

end
