function [c, numSplit] = e_splitdisconnected(c, A)
%E_SPLITDISCONNECTED  Split communities that are not connected in their graph.
%
%   C = PKG.E_SPLITDISCONNECTED(C, A) takes a community index per node and
%   the symmetric adjacency matrix A the communities were found on, and
%   gives each connected piece of a community its own index. Louvain can
%   leave a community in pieces when it moves out the node that bridged
%   them (Traag et al. 2019, DOI:10.1038/s41598-019-41695-z); Leiden
%   guarantees connected communities, and this restores that guarantee
%   after Louvain. Splitting a disconnected community never lowers its
%   modularity, at any resolution.
%
%   When no community is disconnected, C is returned unchanged.
%
%   [C, NUMSPLIT] = PKG.E_SPLITDISCONNECTED(...) also returns how many
%   communities were split.
%
%   See also PKG.E_LOUVAIN, SC_LOUVAIN.

c = c(:);
n = numel(c);
[i, j] = find(A);
keep = c(i) == c(j) & i ~= j;
G = graph(sparse(i(keep), j(keep), 1, n, n), 'upper');
pieces = conncomp(G)';

nc = numel(unique(c));
if max(pieces) == nc
    numSplit = 0;
    return
end

% The pieces never cross communities, so they are the new communities.
pairs = unique([c, pieces], 'rows');
numSplit = sum(accumarray(pairs(:, 1), 1) > 1);
c = pieces;
end
