function [clusterLabels, T, details] = sc_celltypeannoref(X, g, clusterid, ref, options)
%SC_CELLTYPEANNOREF Label clusters against a summary-statistic cell-type reference.
%
%   [clusterLabels, T] = SC_CELLTYPEANNOREF(X, g, clusterid, ref) labels
%   each cluster of the query by comparing it with the cell types of REF,
%   a reference made by SC_BUILDCELLTYPEREF or assembled by hand from
%   published tables. Neither dataset is integrated with the other and no
%   reference cells are needed, so the cost does not grow with the size of
%   the atlas the reference was built from.
%
%   This follows SAHA (Acri et al., bioRxiv 2026,
%   doi:10.64898/2026.07.30.741795), which runs two independent tests:
%
%   CORRELATION (marker-free). Each cluster's mean log-normalised profile
%   and each reference type's profile are z-scored per gene ACROSS the
%   groups of their own dataset, then correlated. Centring within each
%   dataset is what makes this work across platforms: depth, chemistry and
%   the absolute scale of expression cancel, and what is compared is which
%   genes a group expresses more than its neighbours do. It also means the
%   answer depends on which other groups are present -- see Warnings.
%
%   MARKERS (marker-based). Each cluster's positive, specific markers
%   (PKG.I_CLUSTERMARKERS) are tested for over-representation in each
%   type's marker set by a hypergeometric test over the shared genes,
%   Benjamini-Hochberg corrected over every cluster-type pair.
%
%   Inputs:
%     X          genes-by-cells raw counts of the query
%     g          gene names, one per row of X
%     clusterid  cluster of each cell
%     ref        struct with Genes, Types, and Profile and/or Markers;
%                see SC_BUILDCELLTYPEREF
%
%   Name-value options:
%     Method           "both" (default), "correlation" or "markers".
%                      "both" labels by correlation and reports the marker
%                      test beside it, so a label the markers do not
%                      support is visible rather than hidden. A reference
%                      with only one of Profile and Markers uses that one.
%     CorrelationType  "Spearman" (default) or "Pearson"
%     GeneFDR          a shared gene is used for correlation only if its
%                      mean differs between the query clusters by a one-way
%                      ANOVA at this BH-adjusted level (0.01). This is what
%                      keeps noise genes out, and it applies however many
%                      genes there are - a 300-gene panel is filtered the
%                      same way as a 30000-gene transcriptome.
%     NumGenes         at most this many of those, largest F first (2000).
%                      Inf keeps every gene that passes.
%     MinGenes         if fewer genes pass GeneFDR, use this many with the
%                      largest F and warn (20)
%     FineTuneDelta    types scoring within this of a cluster's best
%                      correlation are re-decided on only the genes that
%                      separate them, as SingleR does (0.05). 0 turns it
%                      off. Without it, a type nested in another (its
%                      markers a subset of the other's) is a coin toss.
%                      The Correlation column still reports the chosen
%                      type's correlation, which can be below the best.
%     Alpha            FDR below which a marker match counts (0.05). With
%                      Method="markers", a cluster with no match below it
%                      is "undetermined".
%     NumMarkers, MinLog2FC, MinPctIn, MaxPctOut
%                      query marker filter, passed to PKG.I_CLUSTERMARKERS
%
%   Outputs:
%     clusterLabels  one label per cluster, in the order of T
%     T              one row per cluster: Cluster, NumCells, CellType,
%                    Correlation (of CellType), BestMarkerType, FDR (of
%                    CellType in the marker test), Overlap (share of the
%                    cluster's markers found in CellType's), Agree (the two
%                    tests name the same type)
%     details        Types, Rho (types-by-clusters), FDR (types-by-
%                    clusters), GenesUsed
%
%   Warnings:
%     Fewer than three clusters, or three reference types, make the
%     z-scores degenerate: with two groups every gene scores +-0.71, so the
%     correlation only says which of the two a cluster is closer to. A
%     reference missing the true cell type still returns its nearest type;
%     read Correlation and FDR before trusting a label.
%
%   Example:
%     ref = sc_buildcelltyperef(atlas.X, atlas.g, atlas.c_cell_type_tx);
%     [lab, T] = sc_celltypeannoref(sce.X, sce.g, sce.c_cluster_id, ref);
%     disp(T(~T.Agree, :))            % clusters the two tests disagree on
%
%   See also SC_BUILDCELLTYPEREF, SC_ANNOTATECELLS, SC_SINGLER.

arguments
    X {mustBeNumeric}
    g {mustBeText}
    clusterid
    ref (1, 1) struct
    options.Method (1, 1) string {mustBeMember(options.Method, ...
        ["both", "correlation", "markers"])} = "both"
    options.CorrelationType (1, 1) string {mustBeMember(options.CorrelationType, ...
        ["Spearman", "Pearson"])} = "Spearman"
    options.GeneFDR (1, 1) double {mustBeBetween(options.GeneFDR, 0, 1)} = 0.01
    options.NumGenes (1, 1) double {mustBePositive} = 2000
    options.MinGenes (1, 1) double {mustBePositive, mustBeInteger} = 20
    options.FineTuneDelta (1, 1) double {mustBeNonnegative} = 0.05
    options.Alpha (1, 1) double {mustBeBetween(options.Alpha, 0, 1)} = 0.05
    options.NumMarkers (1, 1) double {mustBePositive} = 100
    options.MinLog2FC (1, 1) double = log2(1.5)
    options.MinPctIn (1, 1) double = 0.5
    options.MaxPctOut (1, 1) double = 0.25
end

g = string(g(:));
if numel(g) ~= size(X, 1)
    error("sc_celltypeannoref:GeneCount", ...
        "g has %d entries for %d rows of X. Pass one gene name per row.", ...
        numel(g), size(X, 1));
end
if numel(clusterid) ~= size(X, 2)
    error("sc_celltypeannoref:ClusterCount", ...
        "clusterid has %d entries for %d cells. Pass one cluster per cell.", ...
        numel(clusterid), size(X, 2));
end
[types, hasProfile, hasMarkers] = in_checkref(ref);
[c, cL] = pkg.i_grp2idxsorted(clusterid);
numClusters = numel(cL);
numTypes = numel(types);

useCorrelation = hasProfile && options.Method ~= "markers";
useMarkers = hasMarkers && options.Method ~= "correlation";
if ~useCorrelation && ~useMarkers
    error("sc_celltypeannoref:MethodNotSupported", ...
        ['Method "%s" needs ref.%s, which is empty. Use Method="both" to ', ...
         'fall back to the test the reference supports.'], options.Method, ...
        in_ifelse(options.Method == "markers", "Markers", "Profile"));
end

rho = nan(numTypes, numClusters);
genesUsed = strings(0, 1);
if useCorrelation
    if numClusters < 3 || numTypes < 3
        warning("sc_celltypeannoref:FewGroups", ...
            ['The correlation test z-scores across groups, and there are ', ...
             '%d clusters and %d reference types. Below three the scores ', ...
             'are degenerate; treat the labels as a two-way split.'], ...
            numClusters, numTypes);
    end
    [rho, genesUsed, Zr, Zq] = in_correlate(X, g, c, numClusters, ref, options);
    bestByRho = in_finetune(rho, Zr, Zq, options.FineTuneDelta);
else
    bestByRho = ones(1, numClusters);
end

fdr = nan(numTypes, numClusters);
overlap = nan(numTypes, numClusters);
if useMarkers
    queryMarkers = pkg.i_clustermarkers(X, g, c, ...
        NumMarkers=options.NumMarkers, MinLog2FC=options.MinLog2FC, ...
        MinPctIn=options.MinPctIn, MaxPctOut=options.MaxPctOut);
    universe = g;
    if hasProfile
        universe = intersect(upper(g), upper(string(ref.Genes(:))));
    end
    [P, ~, overlap] = pkg.e_markeroverlap(queryMarkers, ref.Markers, universe);
    fdr = pkg.e_fdr(P);
end

[minFdr, bestByMarker] = min(fdr, [], 1);
isNoMarkerMatch = isnan(minFdr);
bestMarkerType = types(bestByMarker);
bestMarkerType(isNoMarkerMatch) = "";

if useCorrelation
    chosen = bestByRho;
    clusterLabels = types(chosen);
else
    chosen = bestByMarker;
    clusterLabels = bestMarkerType;
    clusterLabels(isNoMarkerMatch | minFdr >= options.Alpha) = "undetermined";
end

pick = sub2ind([numTypes, numClusters], chosen, 1:numClusters);
isUndetermined = clusterLabels == "undetermined";
corrOfLabel = rho(pick);
fdrOfLabel = fdr(pick);
overlapOfLabel = overlap(pick);
corrOfLabel(isUndetermined) = NaN;
fdrOfLabel(isUndetermined) = NaN;
overlapOfLabel(isUndetermined) = NaN;

agree = useCorrelation & useMarkers & ~isNoMarkerMatch & ...
    (bestByRho == bestByMarker) & (minFdr < options.Alpha);

clusterLabels = clusterLabels(:);
T = table(string(cL(:)), accumarray(c(:), 1, [numClusters, 1]), clusterLabels, ...
    corrOfLabel(:), bestMarkerType(:), fdrOfLabel(:), overlapOfLabel(:), ...
    agree(:), 'VariableNames', ["Cluster", "NumCells", "CellType", ...
    "Correlation", "BestMarkerType", "FDR", "Overlap", "Agree"]);
details = struct("Types", types, "Rho", rho, "FDR", fdr, ...
    "GenesUsed", genesUsed);
end


function [rho, genesUsed, Zr, Zq] = in_correlate(X, g, c, numClusters, ref, options)
[~, iq, ir] = intersect(upper(g), upper(string(ref.Genes(:))), "stable");
if numel(iq) < 10
    error("sc_celltypeannoref:NoSharedGenes", ...
        ['The query and reference share %d genes. Check that both use ', ...
         'the same gene naming (symbols vs Ensembl IDs) and species.'], ...
        numel(iq));
end

Xn = log1p(sc_norm(X));
Xn = Xn(iq, :);
numCells = size(X, 2);
member = sparse(1:numCells, c, 1, numCells, numClusters);
clusterSize = full(sum(member, 1));
queryProfile = full(Xn*member)./clusterSize;
refProfile = double(ref.Profile(ir, :));

% Keep only genes whose cluster means differ by more than within-cluster
% noise explains: a one-way ANOVA across the query clusters, BH-corrected.
% Z-scoring gives every gene equal weight in the correlation, so a gene
% that differs between clusters only by sampling noise adds a random
% +-1 to every profile. Selecting by the variance of the cluster means
% alone, which this once did and only when there were more than NumGenes
% genes, let that noise in whenever the gene list was short: with 80
% marker genes among 880, two types differing in 12 genes were given the
% same label, at correlations of 0.02-0.12.
F = in_anovaf(Xn, queryProfile, clusterSize);
p = fcdf(F, numClusters - 1, numCells - numClusters, "upper");
q = pkg.e_fdr(p);

% A gene constant across the reference's types has no z-score there.
isInformative = q < options.GeneFDR & std(refProfile, 0, 2) > 0;
if nnz(isInformative) < options.MinGenes
    warning("sc_celltypeannoref:FewInformativeGenes", ...
        ['Only %d shared genes differ between the query clusters at FDR ', ...
         '< %g. The correlation uses the %d with the largest F instead; ', ...
         'treat its labels with caution.'], nnz(isInformative), ...
        options.GeneFDR, options.MinGenes);
    F(std(refProfile, 0, 2) == 0) = -Inf;
    [~, order] = sort(F, "descend");
    isInformative = false(size(F));
    isInformative(order(1:min(options.MinGenes, numel(order)))) = true;
end
idx = find(isInformative);
if numel(idx) > options.NumGenes
    [~, order] = sort(F(idx), "descend");
    idx = sort(idx(order(1:options.NumGenes)));
end
queryProfile = queryProfile(idx, :);
refProfile = refProfile(idx, :);
genesUsed = g(iq(idx));

Zr = normalize(refProfile, 2);
Zq = normalize(queryProfile, 2);
rho = corr(Zr, Zq, "type", char(options.CorrelationType));
end


function best = in_finetune(rho, Zr, Zq, delta)
% SingleR-style refinement of near-ties. A type nested in another - its
% markers a subset of the other's - differs from it in a handful of the
% genes the correlation uses, so the two score within noise of each other.
% Among the types within DELTA of a cluster's best correlation, pick the
% one nearest in z-score over only the genes where those candidates differ
% by at least one standard deviation. Distance, not correlation: on those
% genes the nested type is flat, and a flat vector has no correlation.
[~, best] = max(rho, [], 1);
minSeparation = 1;
for k = 1:size(rho, 2)
    cand = find(rho(:, k) >= max(rho(:, k)) - delta);
    if numel(cand) < 2
        continue
    end
    spread = max(Zr(:, cand), [], 2) - min(Zr(:, cand), [], 2);
    genes = spread >= minSeparation;
    if ~any(genes)
        continue
    end
    dist = sum((Zq(genes, k) - Zr(genes, cand)).^2, 1);
    [~, nearest] = min(dist);
    best(k) = cand(nearest);
end
end


function F = in_anovaf(Xn, clusterMeans, clusterSize)
% One-way ANOVA F per gene (row of Xn) across the groups of CLUSTERMEANS,
% computed from per-group sums so that Xn stays sparse.
numCells = size(Xn, 2);
numGroups = numel(clusterSize);
grandMean = full(sum(Xn, 2))/numCells;
ssBetween = ((clusterMeans - grandMean).^2)*clusterSize(:);
ssTotal = full(sum(Xn.^2, 2)) - numCells*grandMean.^2;
ssWithin = max(ssTotal - ssBetween, 0);
F = (ssBetween/(numGroups - 1))./(ssWithin/(numCells - numGroups));
% A gene with no within-group variance: infinitely informative if the
% means differ, uninformative if nothing varies at all.
F(ssWithin == 0 & ssBetween > 0) = Inf;
F(ssBetween == 0) = 0;
end


function [types, hasProfile, hasMarkers] = in_checkref(ref)
if ~isfield(ref, "Types") || isempty(ref.Types)
    error("sc_celltypeannoref:BadReference", ...
        "ref.Types is missing. Build the reference with SC_BUILDCELLTYPEREF.");
end
types = string(ref.Types(:))';
hasProfile = isfield(ref, "Profile") && ~isempty(ref.Profile);
hasMarkers = isfield(ref, "Markers") && ~isempty(ref.Markers);
if hasProfile
    if ~isfield(ref, "Genes") || numel(ref.Genes) ~= size(ref.Profile, 1)
        error("sc_celltypeannoref:BadReference", ...
            "ref.Genes must name every row of ref.Profile.");
    end
    if size(ref.Profile, 2) ~= numel(types)
        error("sc_celltypeannoref:BadReference", ...
            "ref.Profile has %d columns for %d entries of ref.Types.", ...
            size(ref.Profile, 2), numel(types));
    end
end
if hasMarkers && numel(ref.Markers) ~= numel(types)
    error("sc_celltypeannoref:BadReference", ...
        "ref.Markers has %d sets for %d entries of ref.Types.", ...
        numel(ref.Markers), numel(types));
end
if ~hasProfile && ~hasMarkers
    error("sc_celltypeannoref:BadReference", ...
        "The reference has neither a Profile nor Markers to compare against.");
end
end


function out = in_ifelse(condition, a, b)
if condition
    out = a;
else
    out = b;
end
end
