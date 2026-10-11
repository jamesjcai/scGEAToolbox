function [T, P] = ml_alona(X, genelist, clusterid, varargin)
%ML_ALONA Score PanglaoDB cell types against each cluster (CTA scoring).
%
%   T = ML_ALONA(X, GENELIST, CLUSTERID) returns, per cluster, the ranked
%   cell types and their cell-type activity (CTA) scores. Markers come from
%   the HUMAN or MOUSE sheet of assets/PanglaoDB/celltypes.xlsx.
%
%   [T, P] = ML_ALONA(...) also returns P, laid out like T: per cluster the
%   one-sided Fisher p-value and BH-adjusted FDR of each ranked cell type.
%
%   T = ML_ALONA(..., 'species', "mouse") picks the sheet.
%   T = ML_ALONA(..., 'bestonly', false) keeps up to ten rows per cluster
%   instead of the top one.
%   T = ML_ALONA(..., 'alpha', a) sets the FDR above which the winner is
%   reported as "Unknown" (default 0.05). With 'bestonly' false the ranked
%   candidate list is returned as is, for the caller to choose from; the
%   FDRs in P say which of them pass.
%   T = ML_ALONA(..., 'summary', "mean") scores each gene per cluster by
%   its mean instead of the paper's median (default "median"); the Fisher
%   test's expressed call stays median > 0. This is
%   a deviation from the paper, for droplet data with strong ambient RNA:
%   a transcript present in every cell has a high median in every cluster,
%   while a real marker detected in under half a cluster's cells has a
%   median of zero. On example_data/new_example_sce.mat the median calls
%   fibroblast cluster 13 "Enteroendocrine cells" (ambient Gcg, Iapp, Pyy
%   and Sst are in all its cells; 52 of 99 fibroblast markers have median
%   0) and the mean calls it "Fibroblasts".
%
%   PKG.I_CTAMETHODSTEXT writes a Methods paragraph for a run.
%
%   THE METHOD (Franzen, Gan & Bjorkegren 2019, Database baz046, p. 3-4)
%   1. Counts are scaled to library size and log2(x + 1) transformed.
%   2. Each gene is summarised per cluster by its MEDIAN over the cells
%      (or the mean, with 'summary' "mean").
%   3. Marker gene g gets weight w(g) = 1 + sqrt((max f - f(g)) /
%      (max f - min f)), where f(g) is the number of cell types it marks
%      (PKG.E_MARKERWEIGHT), so specific markers count up to twice as much.
%   4. S(j,k) = sum_i Z(k,j,i) * w(i) / N^(1/3) is the CTA score of type j
%      in cluster k, over its N markers. The top score wins.
%   5. A one-sided Fisher's exact (hypergeometric) test asks whether the
%      type's markers are over-represented among the genes the cluster
%      expresses (median > 0). P-values are BH-adjusted; if the winner's
%      adjusted p is above ALPHA it is reported as "Unknown".
%
%   Three details the paper leaves open, and the reading taken here:
%   - N counts the type's markers present in GENELIST. In the paper every
%     annotated gene was in the matrix, so this equals the type's full
%     marker count there; on gene-filtered data a marker that was dropped
%     is not treated as unexpressed.
%   - The Fisher test's universe is GENELIST: marker / not marker by
%     expressed / not expressed.
%   - The BH adjustment is across the cell types tested in one cluster.
%
%   This used to depart from the paper in three ways: the MEAN of
%   ln(x + 1) rather than the median of log2(x + 1), no test, and no
%   "Unknown". On example_data/new_example_sce.mat the switch moved the
%   labels of 3 of 20 clusters (see the commit that made it).
%
%   This used to be two near-identical functions. ML_ALONA read the text
%   files under external/fun_alona_panglaodb and ML_ALONA_NEW the workbook;
%   the scoring loop was the same in both, and on the same input they
%   ranked the types identically, differing only in the fourth significant
%   figure of the score because the two marker tables are not quite the
%   same vintage. The workbook version survived. Two options went with the
%   text-file version, both already unreachable: 'subtype', whose work
%   SC_CSUBTYPEANNO now does against cellsubtypes.xlsx, and zebrafish,
%   which no caller could request -- CLI.CMD_CELLTYPES rejects any species
%   but human and mouse, and the cell-type callbacks all offer those two.
%   The text files stay where they are; PKG.I_GET_PANGLAODBMARKERS,
%   PKG.E_PRIMARYMARKERS, RUN.R_CLUSTERMOLE and others still read them.
%
% https://alona.panglaodb.se/
% https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041
% REF: PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data

if isempty(X) || isempty(genelist)
    k = 1;
    T = table("Unknown", 0, 'VariableNames', ...
        {sprintf('C%d_Cell_Type', k), sprintf('C%d_CTA_Score', k)});
    P = table(NaN, NaN, 'VariableNames', ...
        {sprintf('C%d_P_Value', k), sprintf('C%d_FDR', k)});
    return;
end
if nargin < 3 || isempty(clusterid)
    clusterid = ones(1, size(X, 2));
end
if min(size(clusterid)) ~= 1 || ~isnumeric(clusterid)
    error('CLUSTERID={vector|[]}');
end

p = inputParser;
addRequired(p, 'X', @isnumeric);
addRequired(p, 'genelist', @isstring);
addRequired(p, 'clusterid', @isnumeric);
addOptional(p, 'species', "human", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["human", "mouse"]));
addOptional(p, 'bestonly', true, @islogical);
addParameter(p, 'alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'summary', "median", @(x) (isstring(x) || ischar(x)) && ...
    ismember(lower(string(x)), ["median", "mean"]));
parse(p, X, genelist, clusterid, varargin{:});
species = p.Results.species;
bestonly = p.Results.bestonly;
alpha = p.Results.alpha;
usemean = lower(string(p.Results.summary)) == "mean";

pth = fullfile(cdgea, 'assets', 'PanglaoDB', 'celltypes.xlsx');

if issparse(X)
    try
        X = full(X);
    catch
        disp('Using sparse input--longer running time is expected.');
    end
end
genelist = upper(genelist(:));

Tm = readtable(pth,'FileType','spreadsheet','Sheet',species);
Tw = pkg.e_markerweight(Tm);
wvalu = Tw.Var2;
wgene = string(upper(Tw.Var1));
validG = intersect(genelist, wgene);

% Steps 1-2: normalize on the full matrix (library size has to mean the
% cell's total counts), log2(x + 1), then each gene's median (or mean) per
% cluster.
% All genes, not only markers: the Fisher test needs every gene's
% expressed / not-expressed call.
X = log2(sc_norm(X) + 1);
clusterid = clusterid(:)';
NC = max(clusterid);
% The expressed call for the Fisher test stays on the median either way:
% nearly every gene has a mean above zero, so a mean-based call would make
% every gene "expressed" and the test unable to reject anything.
Zmedian = zeros(numel(genelist), NC);
Z = Zmedian;
for k = 1:NC
    incluster = clusterid == k;
    if any(incluster)
        Zmedian(:, k) = median(X(:, incluster), 2);
        if usemean
            Z(:, k) = mean(X(:, incluster), 2);
        end
    end
end
if ~usemean, Z = Zmedian; end
isexpressed = Zmedian > 0;
numGenes = numel(genelist);
numExpressed = sum(isexpressed, 1);

celltypev = string(Tm.CellType);
markergenev = string(Tm.PositiveMarkers);
numTypes = length(celltypev);
S = zeros(numTypes, NC);
pval = NaN(numTypes, NC);

for j = 1:numTypes
    g = strtrim(strsplit(markergenev(j), ','));
    g = upper(unique(g(strlength(g) > 0)));
    [found, loc] = ismember(g, genelist);
    if ~any(found), continue; end
    g = g(found);
    loc = loc(found);
    [~, wloc] = ismember(g, wgene);
    w = wvalu(wloc);
    numMarkers = numel(loc);

    % Step 4.
    S(j, :) = (w(:)' * Z(loc, :)) ./ nthroot(numMarkers, 3);

    % Step 5: P(at least this many of the type's markers among the
    % expressed genes), drawing numMarkers genes from numGenes of which
    % numExpressed are expressed.
    numHit = sum(isexpressed(loc, :), 1);
    pval(j, :) = hygecdf(numHit - 1, numGenes, numExpressed, numMarkers, 'upper');
end

fdr = NaN(numTypes, NC);
for k = 1:NC
    tested = ~isnan(pval(:, k));
    if any(tested)
        [~, ~, ~, fdr(tested, k)] = pkg.e_fdr_bh(pval(tested, k), alpha, 'pdep', 'no');
    end
end

T = table();
P = table();
for k = 1:NC
    [c, idx] = sort(S(:, k), 'descend');
    names = celltypev(idx);
    pk = pval(idx, k);
    fk = fdr(idx, k);
    if isempty(validG)
        names = "Unknown";
        c = 0;
        pk = NaN;
        fk = NaN;
    elseif bestonly && ~(fk(1) <= alpha)
        names(1) = "Unknown";
    end
    T = [T, table(names, c, 'VariableNames', ...
        {sprintf('C%d_Cell_Type', k), sprintf('C%d_CTA_Score', k)})]; %#ok<AGROW>
    P = [P, table(pk, fk, 'VariableNames', ...
        {sprintf('C%d_P_Value', k), sprintf('C%d_FDR', k)})]; %#ok<AGROW>
end
if size(T, 1) > 10
    T = T(1:10, :);
    P = P(1:10, :);
end
if bestonly && size(T, 1) > 1
    T = T(1, :);
    P = P(1, :);
end
end
