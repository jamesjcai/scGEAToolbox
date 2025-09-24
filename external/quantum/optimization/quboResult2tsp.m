function [order, tourGraph] = quboResult2tsp(result, G)
%QUBORESULT2TSP Convert QUBO Result to TSP form
%
%   ORDER = QUBORESULT2TSP(RESULT) converts the QUBO result to a solution
%   for the TSP problem. The solution is returned as numeric index vector
%   of the order of locations visited.
%
%   [ORDER, TOURGRAPH] = QUBORESULT2TSP(RESULT, G) returns the ORDER as the
%   node names of the specified graph G. Additionally, the optimal tour is
%   returned as a graph.

%   Copyright 2025 The MathWorks, Inc.

arguments
    result (1, 1) {mustBeA(result, "quboResult"), mustBeForTSP}    
    G {mustBeForProblem(result, G)} = []
end

% If there is no best solution, just return empty items
if isempty(result.BestX)
    order = [];
    tourGraph = [];
    return
end

% Find order of cities in route
binx = result.BestX;
N = sqrt(numel(binx));
binx = reshape(binx, N, N);
order = zeros(1, N);
for i = 1:N
    thisIdx = find(binx(i, :)); 
    % binx is a nTimeStep-by-numCities logical array.
    % For a valid tour, each time step should correspond to one city
    if isscalar(thisIdx)
        order(i) = thisIdx;
    else
        break
    end
end

% For a valid tour, each city should only be visited once.
% If a valid tour has not been found then return empty items
if any(~order) || numel(unique(order)) < N
    order = [];
    tourGraph = [];
    return
end

% Create graph of cities in route
if nargin > 1 && ~isempty(G)
    if isempty(G.Nodes)
        G.Nodes = table(string((1:N)'), 'VariableNames', "Name");        
    end
    cities = G.Nodes{:, 1};
    src = cities(order);
    tgt = cities(circshift(order, -1));      
    order = G.Nodes{order, 1};
    idxEdges = findedge(G, src, tgt);
    wts = G.Edges.Weight(idxEdges);
    tourGraph = graph(src, tgt, wts);
else
    tourGraph = [];
end

end

function mustBeGraphOrEmpty(input)

if ~(isempty(input) || isa(input, "graph"))
    error(message("quantum:annealing:quboResult2tsp:mustBeGraphOrEmpty"));
end

end

function mustBeForTSP(result)

numVars = numel(result.BestX);
numQuboCities = sqrt(numVars);
% Number of cities musr be integer, Using the check from mustBeInteger.
% Reproduced here because we don't need the other checks in mustBeInteger
% and we want to throw a custom error message.
isIntegerCities = numQuboCities == fix(numQuboCities);
if ~isIntegerCities
    error(message("quantum:annealing:quboResult2tsp:mustBeForTSP"));
end

end

function mustBeForProblem(result, G)

mustBeGraphOrEmpty(G);
if isempty(G)
    return
end

quantum.internal.optimization.validators.mustBeFullyConnectedGraphWithNoSelfLoops(G);
numVars = numel(result.BestX);
numQuboCities = sqrt(numVars);
if ~isequal(numVars, numnodes(G)^2)
    error(message("quantum:annealing:quboResult2tsp:mustBeForProblem", ...
        numQuboCities, numnodes(G)));
end

end





