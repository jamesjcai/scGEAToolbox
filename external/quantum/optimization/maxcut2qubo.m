function qp = maxcut2qubo(G)
%MAXCUT2QUBO Convert a max-cut problem to QUBO
%
%   QP = MAXCUT2QUBO(G) converts a max-cut problem given by the graph G to
%   an equivalent QUBO formulation.

%   Copyright 2024 The MathWorks, Inc.

% Reference
% [1] A Tutorial on Formulating and Using QUBO Models, Fred Glover and Gary
% Kochenberger and Yu Du, 2019, https://arxiv.org/abs/1811.11538

arguments
    G {mustBeGraph}
end

if any(contains(G.Edges.Properties.VariableNames, "Weight"))
    wts = G.Edges.Weight;
else
    wts = ones(numedges(G), 1);
end
[src, tgt] = findedge(G);

% Disallow self cycles
if any(src == tgt)
    error(message("quantum:annealing:maxcut2qubo:SelfCyclesDisallowed"));
end

% Convert max-cut to QUBO following section 3.2 of [1]
maxNode = max([src(:); tgt(:)]);
Q = full(sparse(src, tgt, wts, maxNode, maxNode));
Q = Q + Q';
c = -sum(Q, 2);
qp = qubo(Q, c);

function mustBeGraph(in)

if ~isa(in, "graph")
    error(message("quantum:annealing:maxcut2qubo:MustBeGraph"));
end












