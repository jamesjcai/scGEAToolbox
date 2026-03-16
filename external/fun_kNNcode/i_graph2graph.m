function [G, T, p] = i_graph2graph(Graph, doplot)

if nargin < 1 || isempty(Graph)
    error('Input Graph is required.');
end
if nargin < 2
    doplot = true;
end

G = graph;
for k = 1:size(Graph, 2)
    g = Graph(:, k);
    for i = 1:length(g)
        nbr = g(i);
        if nbr ~= k && nbr >= 1 && nbr <= size(Graph, 2)
            G = addedge(G, k, nbr);
        end
    end
end

if doplot
    figure;
    p = plot(G);
    [T, ~] = minspantree(G);
    highlight(p, T, 'EdgeColor', 'r')
else
    p = [];
    T = minspantree(G);
end