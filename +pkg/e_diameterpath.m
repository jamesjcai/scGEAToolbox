function [P, s, t] = e_diameterpath(W)

G = graph(W);
d = distances(G);
d(isinf(d)) = 0;
[~, idx] = max(d(:));
[s, t] = ind2sub(size(d), idx);
P = shortestpath(G, s, t);

% allpaths(G,s,t,'MinPathLength',size(W,1), ...
%     'MaxPathLength',size(W,1),'MaxNumPaths',1);
