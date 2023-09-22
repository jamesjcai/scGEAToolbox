function [scnout] = subnetwork(obj, gid)

%assert(all(contains(T.genelist(1:10),scn.g)))
assert(all(contains(gid, obj.g)));
[~, idx] = ismember(gid, obj.g);
SG = subgraph(obj.G, idx);
A = adjacency(SG, 'weighted');
scnout = SingleCellNetwork(A, gid);
