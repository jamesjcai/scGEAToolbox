function i_doublegraphs(G1,G2)
if nargin<2
    G1=WattsStrogatz(500,25,0.15);
    G2=WattsStrogatz(500,25,0.15);
    G1.Nodes.Name = string([1:500]');
    G2.Nodes.Name = string([1:500]');
end
assert(isequal(G1.Nodes.Name,G2.Nodes.Name));
import gui.*
%%
hFig=figure;
h1=subplot(1,2,1);
p1=plot(G1);
h2=subplot(1,2,2);
p2=plot(G2);
p2.XData=p1.XData;
p2.YData=p1.YData;

tb = uitoolbar(hFig);
pt5pickcolr = uipushtool(tb,'Separator','off');				
pt5pickcolr.Tooltip = 'Switch color maps';
pt5pickcolr.ClickedCallback = {@callback_PickColorMap,2};

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
