function [P, overlap, fraction] = e_markeroverlap(queryMarkers, refMarkers, universe)
%E_MARKEROVERLAP Hypergeometric test of marker-set overlap, every pair.
%
%   [P, overlap, fraction] = PKG.E_MARKEROVERLAP(queryMarkers, refMarkers, universe)
%   tests, for every reference set j and query set k, whether the two share
%   more genes than two random sets of the same sizes drawn from UNIVERSE
%   would. Both marker inputs are cell arrays of gene-name vectors; the
%   outputs are numel(refMarkers)-by-numel(queryMarkers).
%
%     P(j, k)        upper-tail hypergeometric p-value, P(X >= overlap)
%     overlap(j, k)  number of shared genes
%     fraction(j, k) overlap as a share of query set k
%
%   UNIVERSE is every gene that could have been a marker on both sides --
%   normally the genes measured in the query, intersected with those of the
%   reference when the reference is a dataset. Genes outside it are dropped
%   from both sets before counting, so an unmeasured reference marker
%   cannot count against a cluster. Gene names are compared ignoring case.
%
%   A pair where either set is empty after that gets P = NaN: no test was
%   possible, which is different from a test that found nothing. Correct P
%   for the number of pairs with PKG.E_FDR, which leaves NaN out of the
%   family.
%
%   This is the marker-based arm of SAHA (Acri et al., bioRxiv 2026,
%   doi:10.64898/2026.07.30.741795).
%
%   See also PKG.I_CLUSTERMARKERS, PKG.E_FDR, SC_CELLTYPEANNOREF.

arguments
    queryMarkers cell
    refMarkers cell
    universe {mustBeText}
end

universe = unique(upper(string(universe(:))));
N = numel(universe);
Q = in_restrict(queryMarkers, universe);
R = in_restrict(refMarkers, universe);

numRef = numel(R);
numQuery = numel(Q);
P = nan(numRef, numQuery);
overlap = zeros(numRef, numQuery);
fraction = nan(numRef, numQuery);
for k = 1:numQuery
    n = numel(Q{k});
    if n == 0
        continue
    end
    for j = 1:numRef
        K = numel(R{j});
        if K == 0
            continue
        end
        x = sum(ismember(Q{k}, R{j}));
        overlap(j, k) = x;
        fraction(j, k) = x/n;
        P(j, k) = hygecdf(x - 1, N, K, n, "upper");
    end
end
end


function S = in_restrict(sets, universe)
S = cell(size(sets));
for i = 1:numel(sets)
    s = unique(upper(string(sets{i}(:))));
    S{i} = s(ismember(s, universe));
end
end
