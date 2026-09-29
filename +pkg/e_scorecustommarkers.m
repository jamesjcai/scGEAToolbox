function [score, types, lists] = e_scorecustommarkers(X, g, c, Tm)
%E_SCORECUSTOMMARKERS Score every cluster against every list of a marker table.
%
%   [score, types, lists] = PKG.E_SCORECUSTOMMARKERS(X, g, c, Tm)
%
%   X      genes-by-cells raw counts
%   g      gene names, one per row of X
%   c      cluster of each cell, integer-coded 1..K
%   Tm     two-column marker table: cell type name, then its genes as one
%          comma separated string (what GUI.I_GETCUSTOMMARKERS returns)
%
%   SCORE is numTypes-by-K. Each marker's mean log-normalised expression per
%   cluster is z-scored across the clusters, and a type's score for a
%   cluster is the weighted mean z of its measured markers, shrunk toward
%   zero by one marker's worth: sum(w.*z)/(sum(w) + mean(w)), which is
%   n/(n+1) times the mean for equal weights. Weights come from
%   PKG.I_MARKERWEIGHTS and favour markers few types share. A type none of
%   whose markers is measured scores -Inf. A positive score means the
%   cluster expresses the list more than the other clusters do; a cluster
%   whose best score is not positive matches no list.
%
%   TYPES and LISTS are the type names and their upper-case marker lists.
%
%   WHY NOT RAW COUNTS. PKG.E_DETERMINECELLTYPE, which this replaces for
%   clusters, sums median raw counts, so one abundant transcript decides
%   everything: on the mouse pancreas example ambient insulin gives every
%   non-beta cell a median of 30-140 insulin counts against 4-53 for its
%   own markers, and it labelled 6 of 10 true cell types "Beta cells".
%   Scoring against the other clusters cancels anything present everywhere.
%   It also means at least two clusters are needed: with one, every z is 0.
%
%   WHY THE SHRINKAGE. A plain mean trusts one marker as much as ten: a
%   macrophage list with two of its three markers missing from the data was
%   left with C1qa alone, whose z of 1.59 in a mixed endocrine cluster beat
%   three concordant Enteroendocrine markers (mean z 1.23) and took 980
%   cells. Excluding such lists instead also threw away the real
%   macrophages, for which C1qa is right. Summing z-scores (Stouffer)
%   rewards list length without limit, so a subtype listed by its few own
%   genes lost its cluster to the parent type listed by many shared ones.
%   n/(n+1) discounts a lone marker by half and a long list hardly at all.
%
%   See also SC_ANNOTATECELLS, PKG.I_MARKERWEIGHTS, GUI.I_GETCUSTOMMARKERS.

arguments
    X {mustBeNumeric}
    g {mustBeText}
    c (:, 1) {mustBeInteger, mustBePositive}
    Tm table
end

[wvalu, wgene, types, markergenev] = pkg.i_markerweights(Tm);
lists = cell(1, numel(markergenev));
for j = 1:numel(markergenev)
    m = strtrim(split(markergenev(j), ","));
    lists{j} = m(strlength(m) > 0);
end

numClusters = max(c);
score = -inf(numel(types), numClusters);
gUpper = upper(string(g(:)));
allMarkers = unique(vertcat(lists{:}));
allMarkers = allMarkers(ismember(allMarkers, gUpper));
if isempty(allMarkers)
    return
end

[~, rowsInX] = ismember(allMarkers, gUpper);
Xn = log1p(sc_norm(X));
Xn = Xn(rowsInX, :);
numCells = size(Xn, 2);
member = sparse(1:numCells, c, 1, numCells, numClusters);
Z = normalize(full(Xn*member)./full(sum(member, 1)), 2);
Z(isnan(Z)) = 0;                                 % constant across clusters

[~, wrow] = ismember(allMarkers, wgene);
w = ones(numel(allMarkers), 1);
w(wrow > 0) = wvalu(wrow(wrow > 0));
for t = 1:numel(types)
    in = ismember(allMarkers, lists{t});
    if any(in)
        score(t, :) = (w(in)'*Z(in, :))/(sum(w(in)) + mean(w(in)));
    end
end
end
