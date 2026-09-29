function [T, bestResolution, bestClusterId, perCell] = sc_annotationstability(sce, options)
%SC_ANNOTATIONSTABILITY Choose a clustering resolution by how its cell-type labels behave.
%
%   [T, bestResolution, bestClusterId, perCell] = SC_ANNOTATIONSTABILITY(sce)
%   clusters the cells at a range of Louvain resolutions, labels every
%   partition, and recommends a resolution in two steps:
%     1. keep the resolutions that show every cell type the sweep finds
%        reliably (type saturation). Coarser than those, two types share a
%        cluster and the minority takes the majority's label.
%     2. of those, take the coarsest whose Agreement is within
%        AgreementTolerance of the best among them. A resolution can show
%        every type and still have one mixed cluster, when each type in it
%        also has a cluster of its own elsewhere; its cells then disagree
%        with their consensus, and step 2 moves past it.
%
%   This follows the resolution sweep of SAHA (Acri et al., bioRxiv 2026,
%   doi:10.64898/2026.07.30.741795), which chooses the resolution by
%   inspecting the annotation across resolutions. It exists because the
%   number of clusters is otherwise a guess: CLUSTERCELLS' own default k is
%   round(NumCells/100, -1), which is 1 for any population under 500 cells.
%
%   WHY NOT THE MOST STABLE RESOLUTION. The obvious criterion - the
%   resolution whose labels best match each cell's consensus across the
%   sweep - is a majority vote, and it fails exactly when it matters. On
%   simulated data (four types, one of them another plus 3-8 genes of its
%   own, eight resolutions from 0.05 to 3), when most resolutions merged
%   the pair the consensus was the merged answer. Over five sizes of that
%   difference, most-agreement picked a resolution labelling 75% of cells
%   correctly in two of them, where 98-99% was available. That is why
%   Agreement only competes among saturated resolutions (step 2): the
%   two-step rule picked the best accuracy available in all five
%   (98-100%).
%
%   WHY NOT SATURATION ALONE. On the mouse pancreas example with a custom
%   marker list, resolution 0.2 already showed all nine types but one of
%   its clusters held 980 cells of five; saturation alone picked it (91.4%
%   correct), the two-step rule picked 0.7 (97.4%; the best of 0.1-1 was
%   97.5%, and it stepped over a dip to 88.8% at 0.5-0.6). With PanglaoDB
%   markers and with a reference built from half the cells, both rules
%   picked 0.7: 99.0% (best 99.1%) and 93.2% (best 93.2%).
%
%   A type counts as found reliably when at least MinCells cells carry it
%   at two or more resolutions, so a label given to one small fragment at
%   one fine resolution cannot drag the recommendation finer.
%
%   SCE is not modified. Apply the answer with
%     sce.c_cluster_id = bestClusterId;
%     sce = sc_annotatecells(sce, Method=..., ...);
%
%   The clustering is Louvain on a shared-nearest-neighbour graph of the
%   principal components, as SCE.CLUSTERCELLS() does by default, including
%   using batch-corrected (Harmony) components when present. The graph is
%   built once and only the resolution varies, so every partition is of the
%   same graph.
%
%   Name-value options:
%     Resolutions      Louvain resolutions to try (0.1:0.1:1)
%     Method           "markers" (PanglaoDB, default), "custommarkers" or
%                      "reference". Only deterministic, cluster-level
%                      methods: an LLM's run-to-run variation would read as
%                      a resolution effect, and the per-cell methods
%                      (SCimilarity, Azimuth) do not depend on clustering.
%     Species          for "markers" ("human")
%     Markers          for "custommarkers": a two-column table of cell type
%                      names and comma separated marker genes
%     Reference        for "reference"; see SC_BUILDCELLTYPEREF
%     ReferenceMethod  for "reference": "both", "correlation" or "markers"
%     Alpha            FDR below which a label counts as marker-supported
%                      (0.05)
%     MinCells         cells a label needs at a resolution to count as a
%                      type found there (10)
%     AgreementTolerance  how far below the best Agreement among saturated
%                      resolutions the pick may be (0.01)
%     Verbose          print progress (true)
%
%   Outputs:
%     T               one row per resolution: Resolution, NumClusters,
%                     NumTypes (reliably found types present, step 1),
%                     Agreement (share of cells whose label is their
%                     consensus label, step 2), Supported (share of cells in
%                     a cluster whose label has FDR < Alpha), Changed
%                     (share of cells whose label differs from the previous
%                     resolution; NaN in row 1), Recommended (logical).
%                     T.Properties.Description holds a paragraph for the
%                     Methods section of a paper, written from the settings
%                     and data of this run (printed too when Verbose).
%     bestResolution  the recommended resolution
%     bestClusterId   the cluster of each cell at that resolution
%     perCell         one row per cell: ConsensusLabel (the label most
%                     resolutions gave it) and Stability (the share that
%                     did). Low stability marks the cells whose type hangs
%                     on the choice of resolution.
%
%   Supported is a weak guide to resolution. The marker test compares each
%   cluster with all other cells, so a type whose markers are all shared
%   with a related type (a subtype, say) has no markers of its own and is
%   never supported, however it is clustered. Read T rather than taking the
%   pick on trust: the recommendation can only be as fine as the finest
%   resolution tried.
%
%   Example:
%     [T, r, cid] = sc_annotationstability(sce, Species="mouse");
%     disp(T)
%     sce.c_cluster_id = cid;
%     sce = sc_annotatecells(sce, Species="mouse");
%
%   See also SC_ANNOTATECELLS, SC_CELLTYPEANNOREF, SC_LOUVAIN.

arguments
    sce (1, 1) SingleCellExperiment
    options.Resolutions (1, :) double {mustBePositive} = 0.1:0.1:1
    options.Method (1, 1) string {mustBeMember(options.Method, ...
        ["markers", "custommarkers", "reference"])} = "markers"
    options.Species (1, 1) string = "human"
    options.Markers table = table()
    options.Reference struct = struct([])
    options.ReferenceMethod (1, 1) string {mustBeMember(options.ReferenceMethod, ...
        ["both", "correlation", "markers"])} = "both"
    options.Alpha (1, 1) double {mustBeBetween(options.Alpha, 0, 1)} = 0.05
    options.MinCells (1, 1) double {mustBePositive} = 10
    options.AgreementTolerance (1, 1) double {mustBeNonnegative} = 0.01
    options.Verbose (1, 1) logical = true
end

if numel(options.Resolutions) < 2
    error('sc_annotationstability:OneResolution', ...
        'Stability needs at least two resolutions to compare. Pass more in Resolutions.');
end
if options.Method == "custommarkers" && (isempty(options.Markers) || ...
        width(options.Markers) < 2)
    error('sc_annotationstability:NoMarkers', ...
        ['Method "custommarkers" needs Markers, a two-column table of cell ', ...
         'type names and comma separated marker genes.']);
end
if options.Method == "reference" && isempty(options.Reference)
    error('sc_annotationstability:NoReference', ...
        'Method "reference" needs Reference, a struct made by SC_BUILDCELLTYPEREF.');
end
resolutions = sort(options.Resolutions);
numRes = numel(resolutions);
numCells = sce.NumCells;

if options.Verbose
    fprintf('[stability] building the neighbour graph...\n');
end
[pcs, usedHarmony] = i_pcs(sce);
[~, info] = sc_louvain(pcs, [], Resolution=resolutions(1));
W = info.Graph;

% A working object that shares X with SCE (copy-on-write) so that labelling
% through SC_ANNOTATECELLS, with its stashing and provenance, never touches
% the caller's object.
work = SingleCellExperiment(sce.X, sce.g);

clusterIds = zeros(numCells, numRes);
labels = strings(numCells, numRes);
supported = false(numCells, numRes);
numClusters = zeros(numRes, 1);
for r = 1:numRes
    c = i_louvainbysize(W, resolutions(r));
    clusterIds(:, r) = c;
    numClusters(r) = max(c);
    if options.Verbose
        fprintf('[stability] resolution %.3g: %d clusters\n', ...
            resolutions(r), numClusters(r));
    end
    if numClusters(r) < 2
        % Nothing to label by cluster; every cell is one group.
        labels(:, r) = "undetermined";
        continue
    end
    work.c_cluster_id = c;
    [work, Ta] = sc_annotatecells(work, Method=options.Method, ...
        Species=options.Species, Reference=options.Reference, ...
        Markers=options.Markers, ...
        ReferenceMethod=options.ReferenceMethod, KeepOld=false, ...
        Verbose=false);
    labels(:, r) = string(work.c_cell_type_tx);
    [~, row] = ismember(string(c), Ta.Cluster);
    supported(:, r) = Ta.FDR(row) < options.Alpha;
end

[consensus, stability] = i_consensus(labels);
agreement = mean(labels == consensus, 1)';
changed = [NaN; mean(labels(:, 2:end) ~= labels(:, 1:end-1), 1)'];
% A type is present at a resolution when at least MinCells cells carry it
% there, and reliable when it is present at two or more resolutions.
% "undetermined" is the absence of a type, not one.
allTypes = unique(labels(labels ~= "undetermined"));
isPresent = false(numel(allTypes), numRes);
for r = 1:numRes
    [~, idx] = ismember(labels(:, r), allTypes);
    counts = accumarray(idx(idx > 0), 1, [numel(allTypes), 1]);
    isPresent(:, r) = counts >= options.MinCells;
end
isReliable = sum(isPresent, 2) >= 2;
numTypes = sum(isPresent(isReliable, :), 1)';

% Among the resolutions that show every reliable type, the coarsest whose
% labels are as settled as any of theirs. Agreement competes only inside
% that set, so a majority of merging resolutions cannot win it; within it,
% it separates a partition that still has a mixed cluster (its cells
% disagree with their consensus) from one that has resolved. Resolutions
% are sorted, so the first match is the coarsest.
isSaturated = numTypes == max(numTypes);
bestAgreement = max(agreement(isSaturated));
best = find(isSaturated & agreement >= bestAgreement - options.AgreementTolerance, 1);
recommended = false(numRes, 1);
recommended(best) = true;

T = table(resolutions(:), numClusters, numTypes, agreement, ...
    mean(supported, 1)', changed, recommended, 'VariableNames', ...
    ["Resolution", "NumClusters", "NumTypes", "Agreement", "Supported", ...
    "Changed", "Recommended"]);
bestResolution = resolutions(best);
bestClusterId = clusterIds(:, best);
perCell = table(consensus, stability, 'VariableNames', ...
    ["ConsensusLabel", "Stability"]);
T.Properties.Description = i_methodstext(options, resolutions, T(best, :), ...
    numCells, size(pcs, 2), usedHarmony, min(2000, numel(sce.g)));

if options.Verbose
    fprintf('[stability] recommended resolution %.3g (%d clusters)\n', ...
        bestResolution, numClusters(best));
    fprintf('\n[stability] Methods text:\n%s\n\n', T.Properties.Description);
end
end


function [pcs, usedHarmony] = i_pcs(sce)
% The components SCE.CLUSTERCELLS clusters on: Harmony's when batch
% correction left one row per cell, otherwise fresh ones from X.
r = sce.struct_cell_reductions;
usedHarmony = isstruct(r) && isfield(r, 'harmony') && ~isempty(r.harmony) && ...
    size(r.harmony, 1) == sce.NumCells;
if usedHarmony
    pcs = r.harmony;
else
    pcs = pkg.e_cellpcs(sce.X, sce.g);
end
end


function c = i_louvainbysize(W, resolution)
% SC_LOUVAIN's numbering, largest cluster first, on a prebuilt graph.
c = pkg.e_louvain(W, resolution);
counts = accumarray(c(:), 1);
[~, order] = sort(counts, 'descend');
relabel = zeros(numel(counts), 1);
relabel(order) = 1:numel(counts);
c = relabel(c(:));
end


function [consensus, stability] = i_consensus(labels)
% The label each cell was given most often, and the share of resolutions
% that gave it. Ties go to the coarser resolution, the earlier column.
numCells = size(labels, 1);
consensus = strings(numCells, 1);
stability = zeros(numCells, 1);
[u, ~, idx] = unique(labels);
idx = reshape(idx, size(labels));
for i = 1:numCells
    counts = accumarray(idx(i, :)', 1, [numel(u), 1]);
    top = max(counts);
    firstTop = find(counts(idx(i, :)) == top, 1);
    consensus(i) = labels(i, firstTop);
    stability(i) = top/size(labels, 2);
end
end


function txt = i_methodstext(options, resolutions, bestRow, numCells, numPCs, ...
        usedHarmony, numHVGs)
% A Methods paragraph for this run. Every number in it is the one used, so
% it stays true when the options change; the steps it names are those of
% PKG.E_CELLPCS, SC_LOUVAIN and the labelling method chosen.
if usedHarmony
    space = sprintf(['the %d Harmony batch-corrected components ' ...
        '(Korsunsky et al., 2019)'], numPCs);
else
    space = sprintf(['%d principal components of the %s most variable ' ...
        'genes (ranked by deviation from the analytic gamma-Poisson ' ...
        'mean-variance curve implied by library size), computed after ' ...
        'library-size normalisation, log1p transformation and per-gene ' ...
        'scaling to unit variance clipped at 10'], numPCs, i_thousands(numHVGs));
end

switch options.Method
    case "markers"
        label = sprintf(['by matching the marker genes of each cluster with ' ...
            'the %s cell-type markers of PanglaoDB (Franzén et al., 2019)'], ...
            options.Species);
    case "custommarkers"
        label = sprintf(['against a user-supplied list of %d cell types and ' ...
            'their marker genes, each marker scored in a cluster against the ' ...
            'other clusters; a cluster expressing no list''s markers above the ' ...
            'others was left undetermined'], height(options.Markers));
    case "reference"
        ref = options.Reference;
        refName = sprintf('a reference of %d cell types', numel(ref.Types));
        if isfield(ref, 'Source') && strlength(string(ref.Source)) > 0
            refName = sprintf('%s (%s)', refName, string(ref.Source));
        end
        correlation = ['Spearman correlation of the mean log-normalised ' ...
            'profile of each cluster with that of each reference type, both ' ...
            'z-scored per gene across the groups of their own dataset, over ' ...
            'the shared genes whose means differed between clusters (one-way ' ...
            'ANOVA, Benjamini-Hochberg adjusted P < 0.01; at most 2,000)'];
        overlap = ['a hypergeometric test of the overlap between the ' ...
            'cluster''s markers and each type''s markers, Benjamini-Hochberg ' ...
            'corrected over all cluster-type pairs'];
        switch options.ReferenceMethod
            case "correlation"
                how = ['by ', correlation];
            case "markers"
                how = ['by ', overlap];
            otherwise
                % "both": correlation labels, the marker test is reported.
                how = ['by ', correlation, ', with ', overlap, ...
                    ' reported beside each label'];
        end
        label = sprintf(['against %s, following SAHA (Acri et al., 2026): ' ...
            'each cluster took the type chosen %s'], refName, how);
    otherwise
        % The arguments block admits only the three methods above.
        label = '';
end

steps = diff(resolutions);
if numel(resolutions) > 2 && all(abs(steps - steps(1)) < 1e-9)
    resText = sprintf('%d resolutions from %g to %g in steps of %g', ...
        numel(resolutions), resolutions(1), resolutions(end), round(steps(1), 6));
else
    resText = sprintf('%d resolutions (%s)', numel(resolutions), ...
        strjoin(compose("%g", resolutions), ", "));
end

versionStr = pkg.i_get_versionnum();
toolbox = 'scGEAToolbox';
if ~isempty(versionStr)
    toolbox = sprintf('scGEAToolbox version %s', versionStr);
end

txt = sprintf(['The clustering resolution was chosen by how the cell-type ' ...
    'annotation behaved across resolutions, following the resolution sweep ' ...
    'of SAHA (Acri et al., 2026). The %s cells were represented by %s. A ' ...
    'shared-nearest-neighbour graph was built on these components (20 ' ...
    'nearest neighbours; edges weighted by the Jaccard index of the two ' ...
    'neighbourhoods and pruned below 1/15), and Louvain community detection ' ...
    '(Blondel et al., 2008) was run on this one graph at %s. At each ' ...
    'resolution the clusters were labelled %s. A cell type was taken as ' ...
    'reliably detected when at least %d cells carried its label at two or ' ...
    'more resolutions. Each cell''s consensus label was the label most ' ...
    'resolutions gave it, and the agreement of a resolution was the share ' ...
    'of cells whose label there matched their consensus. Among the ' ...
    'resolutions at which every reliably detected type was present, the ' ...
    'coarsest was selected whose agreement was within %g of the highest ' ...
    'agreement among them. Agreement was not compared across all ' ...
    'resolutions, because when most resolutions merge two similar types ' ...
    'the consensus is the merged label. The selected resolution was %g, ' ...
    'giving %d clusters and %d cell types (agreement %.1f%%). The analysis ' ...
    'used sc_annotationstability in %s (Cai et al., 2020).'], ...
    i_thousands(numCells), space, resText, label, options.MinCells, ...
    options.AgreementTolerance, bestRow.Resolution, bestRow.NumClusters, ...
    bestRow.NumTypes, 100*bestRow.Agreement, toolbox);
txt = string(txt);
end

function s = i_thousands(n)
% An integer with thousands separators, as a journal prints it.
s = char(string(n));
for k = numel(s)-3:-3:1
    s = [s(1:k), ',', s(k+1:end)];
end
end
