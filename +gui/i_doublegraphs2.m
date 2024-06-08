function [hFig] = i_doublegraphs2(G1, G2, ~)
if nargin < 3, figname = ''; end
if nargin < 2
    G1 = WattsStrogatz(100, 5, 0.15);
    G2 = WattsStrogatz(100, 5, 0.15);
    G1.Nodes.Name = string((1:100)');
    G2.Nodes.Name = string((1:100)');
    G1.Edges.Weight = rand(size(G1.Edges, 1), 1) * 2;
    G2.Edges.Weight = rand(size(G2.Edges, 1), 1) * 2;
end
assert(isequal(G1.Nodes.Name, G2.Nodes.Name));
import gui.*
import ten.*

%%

mfolder = fileparts(mfilename('fullpath'));
load(fullfile(mfolder, ...
    '..', 'resources', 'TFome', 'tfome_tfgenes.mat'), 'tfgenes');

A = blkdiag(adjacency(G1, 'weighted'), adjacency(G2, 'weighted'));
n = matlab.lang.makeUniqueStrings([table2cell(G1.Nodes); table2cell(G2.Nodes)]);
G = digraph(A, n);

% gui.i_singlegraph(G);
hFig = figure;
[~] = drawnetwork(G);
% hFig = gcf;
hFig.Position(3) = hFig.Position(3) * 2.2;

end


function [p] = drawnetwork(G)

A = adjacency(G, 'weighted');
A = ten.e_filtadjc(A, 0.5);
if issymmetric(A)
    G = graph(A, G.Nodes.Name);
else
    G = digraph(A, G.Nodes.Name);
end

w = 1;
p = plot(G);
% layout(p,'force','UseGravity',true)
n = G.numnodes / 2;
p.XData(1:n) = p.XData(n+1:end) - 10;
p.YData(1:n) = p.YData(n+1:end);

%hold on
%plot(p.XData(1:10), p.YData(1:10),'rx')
%plot(p.XData(end-10:end), p.YData(end-10:end),'b+')

cc = repmat([0, 0.4470, 0.7410], G.numedges, 1);
cc(G.Edges.Weight < 0, :) = repmat([0.8500, 0.3250, 0.0980], ...
    sum(G.Edges.Weight < 0), 1);
p.EdgeColor = cc;
p.NodeFontSize = 2 * p.NodeFontSize;
G.Edges.LWidths = abs(w*G.Edges.Weight/max(G.Edges.Weight));
p.LineWidth = G.Edges.LWidths;
end


function h = WattsStrogatz(N, K, beta)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)', 1, K);
t = s + repmat(1:K, N, 1);
t = mod(t-1, N) + 1;

% Rewire the target node of each edge with probability beta
for source = 1:N
    switchEdge = rand(K, 1) < beta;

    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t == source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;

    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s, t);
end
