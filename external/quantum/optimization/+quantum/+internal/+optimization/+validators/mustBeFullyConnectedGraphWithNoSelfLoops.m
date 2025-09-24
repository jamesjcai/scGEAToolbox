function mustBeFullyConnectedGraphWithNoSelfLoops(G)
%MUSTBEFULLYCONNECTEDGRAPHWITHNOSELFLOOPS Validate that value is a fully
%                                         connected graph
%
%   MUSTBEFULLYCONNECTEDGRAPHWITHNOSELFLOOPS(G) throws an error if G is not
%   a fully connected graph without self loops.

%   Copyright 2025 The MathWorks, Inc.

adjMatrix = adjacency(G);
upperAdjMatrixDu = triu(adjMatrix);
numNodes = size(adjMatrix, 1);
fullyConnected = triu(ones(numNodes), 1);
if ~isequal(upperAdjMatrixDu, fullyConnected)
    error(message("quantum:annealing:validators:mustBeFullyConnectedWithNoSelfCycles"));
end

end