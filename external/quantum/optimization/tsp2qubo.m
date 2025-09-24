function qprob = tsp2qubo(D, nvp) 
%TSP2QUBO Convert TSP to QUBO
%
%   QPROB = TSP2QUBO(D) converts the TSP problem defined by distance matrix
%   D into a QUBO. D must be either a symmetric or triangular matrix.
%
%   QPROB = TSP2QUBO(G) converts the TSP problem defined by the fully
%   connected graph G into a QUBO. G must not contain self-loops.
%
%   QPROB = TSP2QUBO(___, Penalty=VALUE) specifies the scalar penalty
%   factor for the TSP constraints that are incorporated into the
%   QUBO objective.
% 
%   The TSP is converted to QUBO following section 4.3 of [1].
%
%   Reference:
%   [1] A Hybrid Solution Method for the Capacitated Vehicle Routing
%   Problem Using a Quantum Annealer, Feld, S. et al, Frontiers in ICT,6:13
%   (2019). https://arxiv.org/pdf/1811.07403.pdf.

%   Copyright 2023-2025 The MathWorks, Inc.

arguments
    D {mustBeDistanceMatrixOrGraph}
    nvp.Penalty {quantum.internal.optimization.validators.mustBeDouble, mustBeReal, mustBeScalarOrEmpty, mustBeNonempty, mustBeFinite} = 1
end

if isa(D, "graph")
    D = distances(D);
end

if ~issymmetric(D)
     D = (D + D.')/2;
end

% Calculate objective
HB = calculateObjective(D);

% Calculate penalized constraints 
[HA, const] = calculatePenalizedConstraints(D, nvp.Penalty);

% Form penalized objective
Q = HA + HB;

% Create qubo
qprob = qubo(Q, [], const);

end

function HB = calculateObjective(D)
% This function creates a quadratic coefficent matrix that corresponds to
% the objective defined by equation (12) below.

% The problem variables, x(i,j), are 1 if node i is visited at time j and 0
% otherwise. Following section 4.3 of [1], equation (12) defines the
% objective HB as
%
%                          n 
% HB =     sum     D_ui * sum  x_uj * x_i,j+1     (12)
%       (ui in E)        j = 1
%
% where E is the graph that corresponds to all the possible connections
% in the TSP. 
% 
% Note that the wrap around term (i.e return to start city) means that 
% x_i,n+1 = x_i,1

% See the "Convert to QUBO: Theory" and "Convert to QUBO: Code" sections in
% the following documentation example for more explanation
% https://mathworks.com/help/matlab/math/quantum-tsp.html
n = size(D, 1);
B = diag(ones(n-1,1),1);
B(1,n) = B(1,n) + 1; % wrap around term
HB = kron(D,B);
HB = (HB + HB')/2;

end


function [HA, offset] = calculatePenalizedConstraints(D, penalty)
% This function creates a quadratic coefficent matrix that corresponds to
% the penalized constraints defined by equation (11) below.

% The problem variables, x(i,j), are 1 if node i is visited at time j and 0
% otherwise. Following section 4.3 of [1], equation (11) defines the
% penalized constraints HA as
%
%           n         n                   n         n
% HA = A * sum  (1 - sum  x_ij)^2  + A * sum  (1 - sum  x_ij)^2    (11)
%          i=1       j=1                 j=1       i=1 
 
% See the "Convert to QUBO: Theory" and "Convert to QUBO: Code" sections in
% the following documentation example for more explanation
% https://mathworks.com/help/matlab/math/quantum-tsp.html

% Number of nodes
n = size(D, 1);

% First set of constraints (each node visted once)
Aeq1 = kron(eye(n), ones(1,n));
beq1 = ones(n,1);

% Second set of constraints (each time has one node)
Aeq2 = zeros(n,n^2);
for i = 1:n
    Aeq2(i,i:n:n^2) = 1;
end
beq2 = ones(n,1);

% Concatenate constraints
Aeq = [Aeq1; Aeq2];
beq = [beq1; beq2];

% Ensure there is a penalty value per constraint
numAeq = size(Aeq, 1);
penalty = max(D(:))*penalty;
if isscalar(penalty)
    penalty= penalty*ones(numAeq, 1);
end

% Calcuate penalized constraints
HA = sparse(n^2, n^2);
for i = 1:numAeq
    % Quadratic elements 
    HA = HA + penalty(i)*kron(Aeq(i, :),Aeq(i, :)');   
end

% Linear elements
HA = HA + diag(sum(-2*Aeq.*beq.*penalty));

% Constant offset ignored by QUBO representation
offset = sum(penalty.*(beq.^2));

end


function mustBeDistanceMatrixOrGraph(input)

if isa(input, "graph")
    quantum.internal.optimization.validators.mustBeFullyConnectedGraphWithNoSelfLoops(input);
elseif isa(input, "digraph")
    error(message("quantum:annealing:tsp2qubo:mustBeGraph"));
else
    mustBeDistanceMatrix(input);
end

end

function mustBeSquare(Q)
numOfRows = size(Q,1);
numOfCols = size(Q,2);
if numOfRows ~= numOfCols
    error(message("quantum:annealing:tsp2qubo:mustBeSquareMatrix"));
end
end

function mustBeSymmetricOrTriangular(Q)
if ~(issymmetric(Q) || istril(Q) || istriu(Q))
    error(message("quantum:annealing:tsp2qubo:mustBeSymmetricOrTriangular"));
end
end

function mustBeZeroDiagonal(Q)

if any(diag(Q))
    error(message("quantum:annealing:tsp2qubo:mustBeZeroDiagonal"));
end

end

function mustBeDistanceMatrix(Q)

quantum.internal.optimization.validators.mustBeDouble(Q);
mustBeReal(Q);
mustBeSquare(Q);
mustBeNonempty(Q);
mustBeFinite(Q);
mustBeNonnegative(Q);
mustBeSymmetricOrTriangular(Q);
mustBeZeroDiagonal(Q);

end

