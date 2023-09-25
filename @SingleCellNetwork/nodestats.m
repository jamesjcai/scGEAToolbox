function [T, X] = nodestats(obj)

import ten.*
A0 = obj.A;
A0sym = 0.5 * (A0 + A0');
A0sym = ten.e_filtadjc(A0sym);
A0 = ten.e_filtadjc(A0);
genelist=obj.g;
G0 = digraph(A0, genelist, 'omitselfloops');
G0x = graph(A0sym, genelist, 'omitselfloops');

dig_indegree = G0.centrality('indegree');
dig_outdegree = G0.centrality('outdegree');
dig_betweenness = G0.centrality('betweenness');
dig_pagerank = G0.centrality('pagerank');
dig_incloseness = G0.centrality('incloseness');
dig_outcloseness = G0.centrality('outcloseness');
dig_authorities = G0.centrality('authorities');
sym_degree = G0x.centrality('degree');
sym_closeness = G0x.centrality('closeness');
sym_betweeness = G0x.centrality('betweenness');
sym_pagerank = G0x.centrality('pagerank');
sym_eigenvector = G0x.centrality('eigenvector');

genes = genelist;

%%
T = table(genes, dig_indegree, dig_outdegree, dig_betweenness, ...
    dig_pagerank, dig_incloseness, dig_outcloseness, ...
    dig_authorities, sym_degree, sym_closeness, sym_betweeness, ...
    sym_pagerank, sym_eigenvector);

X = table2array(T(:, 2:end));
end
