function [T, info] = enrich(X, genelist, grp, opts)
%GLY.ENRICH  Enrichment/depletion of glyco-module expressing cells.
%
%   T = GLY.ENRICH(X, genelist, grp) asks, for every glyco module and
%   every group of cells, whether that group holds more or fewer cells
%   EXPRESSING the module than the rest of the data does. A cell counts as
%   expressing a module when at least one (see MinDetected) of the module's
%   genes is detected in it. The answer is a log2 odds ratio from a 2-by-2
%   contingency table, with a Fisher exact test and a Benjamini-Hochberg
%   adjusted p-value.
%
%   USAGE:
%     T = gly.enrich(sce.X, sce.g, sce.c_cell_type_tx);
%     T = gly.enrich(X, genelist, grp, MinDetected=2);
%     [T, info] = gly.enrich(X, genelist, grp);   % matrices for plotting
%
%   WHY THIS EXISTS, AND WHEN TO PREFER IT TO GLY.STATE.
%
%   GLY.STATE scores a cell against a module with AddModuleScore, which
%   subtracts the mean of a control gene set drawn from the SAME dataset and
%   matched on expression bin. That makes it a within-dataset contrast: it
%   is built to rank cells inside one experiment, and its population mean is
%   pinned near zero by construction. Measured on GSE218310, the mean score
%   of every module in each of three runs sat within +/-0.15 of zero, which
%   made a between-population comparison of module scores uninformative -
%   the quantity simply does not carry a level that survives leaving the
%   dataset.
%
%   The fraction of expressing cells does. It is an absolute proportion, it
%   needs no control set, and it is the metric recommended for lowly
%   expressed genes, which glycogenes are (Chrysinas et al., NAR Genomics
%   and Bioinformatics 2024, 6:lqae169, who applied it across Tabula
%   Sapiens; see also their note that detection fraction beats mean
%   expression for this gene class). Use GLY.ENRICH when the comparison
%   crosses datasets, batches or experiments; use GLY.STATE when ranking
%   cells within one.
%
%   The two are complementary rather than competing, and they answer
%   different questions: GLY.STATE asks how strongly a cell runs a
%   module, GLY.ENRICH asks how many cells run it at all. A module can
%   move on one and not the other - a pathway switched on in more cells at
%   the same per-cell level shows up here and not there.
%
%   THE ONE THING THAT WILL FOOL YOU: SEQUENCING DEPTH. Detection is a
%   function of how deeply a cell was sequenced, so this metric is only
%   comparable between groups that were sequenced comparably. Measured on
%   GSE218310: the two cultured lines sit at ~21 000 median UMI and ~110
%   detected glycogenes per cell, the PBMC run at 4612 and 34. Every module
%   therefore reads as "depleted" in PBMC, and most of that gap is the
%   depth, not the biology. INFO.medianGenesDetected reports the per-group
%   depth and a warning fires when the groups differ by more than
%   DepthWarnRatio, but the check is a tripwire and not a correction:
%   nothing here rescales for depth. If the groups are unbalanced, subsample
%   the counts to a common depth before calling, or restrict the comparison
%   to groups within one run. The same caution applies between cell types
%   inside one experiment, since RNA content differs by cell type too.
%
%   INPUTS:
%     X          - genes-by-cells expression matrix. Only whether an entry
%                  is nonzero is used, so raw counts, normalized counts and
%                  log-normalized counts all give the same answer. Do not
%                  pass imputed or denoised values: imputation invents
%                  detection, which is the whole measurement here.
%     genelist   - G-by-1 gene symbols, length = rows of X
%     grp        - N-by-1 group label per cell (string, cellstr, categorical
%                  or numeric); cell type, tissue, condition, cluster
%     opts.MinDetected (1)   genes of a module that must be detected in a
%                  cell for it to count as expressing. 1 reproduces the
%                  published definition (mean pathway expression nonzero).
%                  Raise it for large modules, where 1 makes almost every
%                  cell an expressing cell and the odds ratio goes flat.
%     opts.MinGenes    (3)   modules with fewer than this many genes present
%                  in the data are dropped, as in GLY.STATE
%     opts.GeneSets    ([])  table of gene sets with columns Name and Genes
%                  (comma-separated), the shape GLY.GENESETS returns.
%                  Default is the curated glycobiology collection.
%     opts.DepthWarnRatio (1.5)  warn when the most and least deeply
%                  sequenced groups differ by more than this ratio in median
%                  detected genes. Inf silences the check.
%
%   OUTPUT T: one row per (module, group), sorted by q then by |log2OR|
%     module, group        - what was tested
%     nExpressing, nGroup  - expressing cells in the group, cells in it
%     fracInGroup          - nExpressing / nGroup
%     fracOutGroup         - the same fraction over every other cell
%     log2OR               - log2 of (a*d)/(b*c); +/-Inf when a margin is 0
%     p, q                 - Fisher exact, and BH over all rows
%     a, b, c, d           - the contingency table itself, so a degenerate
%                            odds ratio can be read rather than guessed at
%
%   OUTPUT info: .modules, .groups, .log2OR, .q, .fracExpressing (modules by
%     groups matrices for heatmaps), .fracOverall (per module, over all
%     cells - the dataset-independent number), .ncommon (genes found),
%     .medianGenesDetected (per group), and .rhoDepth - per module, the
%     Spearman correlation of its expressing fraction with group depth.
%     READ .rhoDepth BEFORE BELIEVING ANY ENRICHMENT. It needs three or
%     more groups; below that it is NaN.
%
% see also: GLY.STATE, GLY.GENESETS, SC_DEGMAST, PKG.E_FDR

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    grp (:, 1)
    opts.MinDetected (1, 1) double {mustBePositive} = 1
    opts.MinGenes (1, 1) double = 3
    opts.GeneSets = []
    opts.DepthWarnRatio (1, 1) double {mustBePositive} = 1.5
end

if numel(genelist) ~= size(X, 1)
    error("GLY:ENRICH:GeneCount", ...
        "GENELIST length (%d) must equal the number of rows of X (%d).", ...
        numel(genelist), size(X, 1));
end
nCells = size(X, 2);
if numel(grp) ~= nCells
    error("GLY:ENRICH:GroupCount", ...
        "GRP length (%d) must equal the number of columns of X (%d).", ...
        numel(grp), nCells);
end

[setnames, setgenelists] = i_getsets(opts.GeneSets);
nSets = numel(setnames);

% Detection per module per cell. Sparsity and single precision both have to
% go: SCE.X is single sparse from R2025a on, and a sparse logical index is
% not what the contingency counts below want.
expressing = false(nSets, nCells);
ncommonAll = zeros(nSets, 1);
upperlist = upper(genelist);
for k = 1:nSets
    idx = ismember(upperlist, upper(setgenelists{k}));
    ncommonAll(k) = sum(idx);
    if ncommonAll(k) < opts.MinGenes
        continue;
    end
    ndet = full(double(sum(X(idx, :) > 0, 1)));
    expressing(k, :) = ndet >= opts.MinDetected;
end

keep = ncommonAll >= opts.MinGenes;
expressing = expressing(keep, :);
modules = setnames(keep);
ncommon = ncommonAll(keep);
nSets = numel(modules);
if nSets == 0
    error("GLY:ENRICH:NoModules", ...
        "No module had at least %d genes present in the data.", opts.MinGenes);
end

[gidx, groups] = findgroups(string(grp(:)));
nGroups = numel(groups);
if nGroups < 2
    error("GLY:ENRICH:OneGroup", ...
        "GRP has a single level; there is nothing to compare it against.");
end

% Depth tripwire. Detection rises with sequencing depth, so an odds ratio
% between groups of unequal depth is partly a statement about the library
% and not about glycosylation. Computed on ALL genes, not just the modules,
% because that is the quantity that is confounded.
genesPerCell = full(double(sum(X > 0, 1)))';
medDepth = accumarray(gidx, genesPerCell, [nGroups, 1], @median);
if isfinite(opts.DepthWarnRatio) && min(medDepth) > 0
    ratio = max(medDepth) / min(medDepth);
    if ratio > opts.DepthWarnRatio
        [~, hi] = max(medDepth);
        [~, lo] = min(medDepth);
        warning("GLY:ENRICH:DepthImbalance", ...
            ['Groups differ %.1f-fold in sequencing depth ("%s" %.0f genes ', ...
            'per cell, "%s" %.0f). Detection scales with depth, so these ', ...
            'odds ratios are confounded: the shallower group will look ', ...
            'depleted in every module. Subsample to a common depth, or ', ...
            'compare only within a run.'], ...
            ratio, groups(hi), medDepth(hi), groups(lo), medDepth(lo));
    end
end

log2OR = NaN(nSets, nGroups);
pval = NaN(nSets, nGroups);
fracIn = NaN(nSets, nGroups);
fracOut = NaN(nSets, nGroups);
A = zeros(nSets, nGroups); B = A; C = A; D = A;

nInGroup = accumarray(gidx, 1, [nGroups, 1])';
for k = 1:nSets
    e = expressing(k, :);
    % Expressing cells per group, and the totals, without rebuilding masks.
    expPerGroup = accumarray(gidx(e), 1, [nGroups, 1])';
    totExp = sum(e);
    for j = 1:nGroups
        a = expPerGroup(j);                 % expressing, in group
        b = nInGroup(j) - a;                % not expressing, in group
        c = totExp - a;                     % expressing, outside
        d = (nCells - nInGroup(j)) - c;     % not expressing, outside
        A(k, j) = a; B(k, j) = b; C(k, j) = c; D(k, j) = d;
        fracIn(k, j) = a / max(nInGroup(j), 1);
        fracOut(k, j) = c / max(nCells - nInGroup(j), 1);
        % log2 of (a*d)/(b*c). A zero margin is reported as +/-Inf rather
        % than smoothed away: the counts are in the table, and a module
        % expressed in every cell of every group is a fact about the module,
        % not something a continuity correction should hide.
        log2OR(k, j) = log2(a) + log2(d) - log2(b) - log2(c);
        pval(k, j) = i_fisher([a, b; c, d]);
    end
end

q = pkg.e_fdr(pval);

[mi, gi] = ndgrid(1:nSets, 1:nGroups);
T = table(modules(mi(:)), groups(gi(:)), A(:), nInGroup(gi(:))', ...
    fracIn(:), fracOut(:), log2OR(:), pval(:), q(:), ...
    A(:), B(:), C(:), D(:), ...
    VariableNames = ["module", "group", "nExpressing", "nGroup", ...
    "fracInGroup", "fracOutGroup", "log2OR", "p", "q", "a", "b", "c", "d"]);
[~, ord] = sortrows([q(:), -abs(log2OR(:))], [1, 2], "ascend", ...
    MissingPlacement = "last");
T = T(ord, :);

info = struct();
info.modules = modules;
info.groups = groups;
info.log2OR = log2OR;
info.p = pval;
info.q = q;
info.fracExpressing = fracIn;
info.fracOverall = sum(expressing, 2) / nCells;
info.ncommon = ncommon;
info.nInGroup = nInGroup';
info.minDetected = opts.MinDetected;
info.medianGenesDetected = medDepth;

% Per module, how strongly its expressing fraction tracks group depth. This
% is the number that separates a real enrichment from a depth artifact, and
% it is not optional reading: on GSE218310 PBMC across 80 clusters the
% median module scored 0.80 here, meaning most of what looks like biology is
% the library. The modules worth trusting are the ones that stand out as
% LOW. Glyco_siglecs at 0.53 was 25th of 26, and separately the clusters at
% exactly 0% expressing turned out no shallower than the rest - which is
% what confirmed it rather than this correlation alone.
info.rhoDepth = NaN(nSets, 1);
if nGroups >= 3
    for k = 1:nSets
        info.rhoDepth(k) = corr(fracIn(k, :)', medDepth, ...
            Type = "Spearman", Rows = "pairwise");
    end
end

end


% ----------------------------------------------------------------------
function [names, genelists] = i_getsets(G)
if isempty(G)
    [setmatrx, names, setgenes] = gly.genesets();
    genelists = cell(numel(names), 1);
    for k = 1:numel(names)
        genelists{k} = setgenes(setmatrx(k, :));
    end
    return;
end
if ~istable(G) || ~all(ismember(["Name", "Genes"], G.Properties.VariableNames))
    error("GLY:ENRICH:BadGeneSets", ...
        "GeneSets must be a table with Name and Genes columns.");
end
names = string(G.Name);
genelists = cell(numel(names), 1);
for k = 1:numel(names)
    v = strtrim(split(string(G.Genes(k)), ","));
    genelists{k} = v(strlength(v) > 0);
end
end


function p = i_fisher(tbl)
% Fisher exact where it is affordable, chi-square where it is not.
%
% FISHERTEST enumerates the hypergeometric tail, and its cost grows with the
% smallest margin. On an atlas - hundreds of thousands of cells, a module
% expressed in most of them - that margin is large enough for the exact test
% to dominate the whole run, while the chi-square approximation it would be
% agreeing with to many digits costs nothing. The switch is on the smallest
% margin rather than on N, because that is what actually drives the cost and
% also what governs when the approximation is safe.
if any(tbl(:) < 0) || any(isnan(tbl(:)))
    p = NaN;
    return;
end
n = sum(tbl(:));
if n == 0
    p = NaN;
    return;
end
rowsum = sum(tbl, 2);
colsum = sum(tbl, 1);
if any(rowsum == 0) || any(colsum == 0)
    % A degenerate table carries no evidence of association.
    p = 1;
    return;
end
if min([rowsum(:); colsum(:)]) <= 1e4
    % FISHERTEST returns the hypothesis DECISION first and the p-value
    % second. Taking the first output silently yields a 0/1 logical that
    % looks like a p-value and passes straight through BH.
    [~, p] = fishertest(tbl);
    return;
end
expected = rowsum * colsum / n;
% Yates-corrected chi-square on 1 df.
chi2 = sum((abs(tbl - expected) - 0.5).^2 ./ expected, "all");
p = gammainc(chi2 / 2, 0.5, "upper");
end
