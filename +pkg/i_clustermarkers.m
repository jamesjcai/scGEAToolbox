function markers = i_clustermarkers(X, g, c, options)
%I_CLUSTERMARKERS Positive, specific marker genes of every group, as gene sets.
%
%   markers = PKG.I_CLUSTERMARKERS(X, g, c) returns a 1-by-K cell array
%   whose k-th entry is a string column of the marker genes of group k,
%   where C is integer-coded 1..K. Each group is tested against all other
%   cells by PKG.E_FINDALLMARKERS (Wilcoxon on log1p-normalised counts),
%   then filtered for specificity and ranked.
%
%   The defaults follow the marker filter of SAHA (Acri et al., bioRxiv
%   2026, doi:10.64898/2026.07.30.741795): adjusted p < 0.05, fold change
%   > 1.5, detected in fewer than 25% of the other cells, top 100 per
%   group. They are strict on purpose: a marker set is used here for an
%   overlap test, and a set padded with weak, broadly expressed genes
%   overlaps with everything.
%
%   The one departure is MinPctIn, 0.5 where SAHA uses 0.75. Detection
%   rate falls with sequencing depth, so 0.75 empties the sets of shallow
%   data. Measured on example_data/new_example_sce.mat (mouse pancreas,
%   10 types), reference from one random half, the other half as query:
%   at 0.75 the Enteroendocrine type got no markers at all, and the marker
%   test named 8/10 types at full depth and 5/10 with the query thinned to
%   25% of its counts. At 0.5 it named 10/10 and 9/10.
%
%   Name-value options:
%     NumMarkers   genes kept per group, strongest first      (100)
%     MinLog2FC    minimum avg_log2FC                         (log2(1.5))
%     MinPctIn     minimum detection rate inside the group    (0.5)
%     MaxPctOut    maximum detection rate outside the group   (0.25)
%     Alpha        cut-off on the Bonferroni-adjusted p-value (0.05)
%
%   A group with no gene passing every filter gets an empty set; the
%   overlap test then reports NaN for it rather than a p-value.
%
%   See also PKG.E_FINDALLMARKERS, PKG.E_MARKEROVERLAP.

arguments
    X {mustBeNumeric}
    g {mustBeText}
    c (:, 1) {mustBeInteger, mustBePositive}
    options.NumMarkers (1, 1) double {mustBePositive} = 100
    options.MinLog2FC (1, 1) double = log2(1.5)
    options.MinPctIn (1, 1) double {mustBeBetween(options.MinPctIn, 0, 1)} = 0.5
    options.MaxPctOut (1, 1) double {mustBeBetween(options.MaxPctOut, 0, 1)} = 0.25
    options.Alpha (1, 1) double {mustBeBetween(options.Alpha, 0, 1)} = 0.05
end

numGroups = max(c);
cL = string(1:numGroups)';
T = pkg.e_findallmarkers(X, string(g(:)), c, cL, options.MinLog2FC, [], ...
    false, [], OnlyPos=true, ReturnThresh=options.Alpha, ...
    ThresholdOn="p_val_adj");
T = T(T.pct_1 >= options.MinPctIn & T.pct_2 <= options.MaxPctOut, :);

markers = cell(1, numGroups);
for k = 1:numGroups
    t = T(T.grp == cL(k), :);
    t = sortrows(t, ["p_val_adj", "avg_log2FC"], ["ascend", "descend"]);
    markers{k} = string(t.g(1:min(height(t), options.NumMarkers)));
end
end
