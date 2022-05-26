function [hFig]=i_doublegraphs2(G1,G2,figname)
if nargin<3, figname=''; end
if nargin<2
    G1=WattsStrogatz(100,5,0.15);
    G2=WattsStrogatz(100,5,0.15);
    G1.Nodes.Name = string((1:100)');
    G2.Nodes.Name = string((1:100)');
    G1.Edges.Weight=rand(size(G1.Edges,1),1)*2;
    G2.Edges.Weight=rand(size(G2.Edges,1),1)*2;
end
assert(isequal(G1.Nodes.Name,G2.Nodes.Name));
import gui.*
import ten.*
%%

mfolder=fileparts(mfilename('fullpath'));
load(fullfile(mfolder,...
       '../resources','tfome_tfgenes.mat'),'tfgenes');

A=blkdiag(adjacency(G1),adjacency(G2));
n=matlab.lang.makeUniqueStrings([table2cell(G1.Nodes); table2cell(G2.Nodes)]);
G=digraph(A,n);
gui.i_singlegraph(G);

end




function h = WattsStrogatz(N,K,beta)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
end
