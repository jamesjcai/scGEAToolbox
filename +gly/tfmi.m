function [M, info] = tfmi(X, genelist, opts)
%GLY.TFMI  Mutual information between TFs and glycogenes.
%
%   [M, info] = GLY.TFMI(X, genelist) scores every transcription factor
%   against every glycogene by the mutual information of their single-cell
%   expression, and returns the TFs-by-glycogenes matrix M. High MI says the
%   two move together in a way that need not be monotone, which is the case
%   for TF-target relations with thresholds or saturation in them.
%
%   USAGE:
%     [M, info] = gly.tfmi(sce.X, sce.g);
%     M = gly.tfmi(sce.X, sce.g, Species="mm", NumBins=4);
%     info.topTFs(1:10, :)      % TFs ranked by total MI over glycogenes
%
%   Literature reports few experimentally validated TFs of glycogenes, and
%   this is the gap Chrysinas et al. (NAR Genomics and Bioinformatics 2024,
%   6:lqae169) set out to fill across Tabula Sapiens. They report MI
%   recovering TFLink's TF-target edges at AUPRC 0.447 against 0.396 for
%   Spearman, 0.367 for Pearson and 0.318 for a random predictor.
%
%   VALIDATE THE RUN BEFORE READING IT. That result is theirs, on 456101
%   cells. Whether MI beats a correlation on YOUR data is a question about
%   your cell count and depth, not a settled fact, and this function
%   answers it rather than assuming it: with Validate=true (the default)
%   info.auprc scores M, |Pearson| and |Spearman| against the DoRothEA
%   TF-target network, alongside the prevalence a random predictor would
%   reach.
%
%   MEASURED HERE, MI LOSES. On GSE218310, against DoRothEA, at every cell
%   count available:
%
%     PBMC       1000 cells  MI 0.0213  Pearson 0.0232  Spearman 0.0228
%     PBMC       3000        MI 0.0228  Pearson 0.0238  Spearman 0.0239
%     PBMC       8344        MI 0.0232  Pearson 0.0239  Spearman 0.0240
%     Panc04.03  5838        MI 0.0211  Pearson 0.0222  Spearman 0.0213
%     (prevalence 0.0202-0.0215)
%
%   MI does climb with cell count, as an entropy estimator should, but it
%   never catches either correlation, and the deeper of the two datasets
%   does not rescue it. Compare lift over prevalence rather than raw AUPRC
%   across studies - the paper's grid against TFLink had prevalence 0.318
%   and ours against DoRothEA has 0.021, so the absolute numbers are not
%   comparable. Their best is 1.40x prevalence; the best here is 1.12x.
%
%   AND WHAT IT DOES FIND IS CELL TYPE, NOT REGULATION. The strongest pairs
%   on the PBMC run are SPI1-NAGK, SPI1-HK3, FOS-NAGK. SPI1 is detected in
%   0-100% of cells depending on cluster, HK3 in 0-76%, NAGK in 4-97%: they
%   are myeloid identity genes, and in a heterogeneous population any two
%   markers of one cell type co-vary strongly whether or not either
%   regulates the other. On this evidence the matrix is a co-expression
%   readout dominated by composition. It is worth computing, and it is not
%   worth clustering into regulatory modules - a matrix that encodes cell
%   identity will cluster into confident-looking modules that are cell
%   types wearing a regulatory label. Reach for SC_TFACTIVITY, which uses
%   DoRothEA's edges rather than trying to rediscover them, unless you have
%   the cell numbers to do better.
%
%   DISCRETIZATION. MI needs binned data, and scRNA-seq counts are
%   zero-inflated, so equal-frequency binning collapses onto the zeros.
%   Bin 1 is therefore "not detected" and the nonzero values are split into
%   NumBins-1 quantile bins. A gene detected in too few cells cannot support
%   an estimate at all and is dropped; see MinDetected.
%
%   INPUTS:
%     X        - genes-by-cells counts, raw or normalized
%     genelist - G-by-1 gene symbols
%     opts.Species     ("hs")  TF list, via PKG.E_GETTFLIST
%     opts.GeneSets    ([])    glycogene collection; default
%                      GLY.ENZONTO, all terms
%     opts.NumBins     (5)     bins, including the zero bin
%     opts.MaxCells    (3000)  subsample; MI is quadratic in the bin grid
%                      and linear in cells, and 3000 is ample for 5 bins
%     opts.MinDetected (0.02)  drop genes and TFs detected in a smaller
%                      fraction of cells than this
%     opts.Validate    (true)  compute info.auprc
%     opts.Seed        (42)
%
%   OUTPUT info: .tfs, .glycogenes, .nCells, .auprc (a table of estimator
%     against AUPRC when Validate), .topTFs, .topPairs, .fracDetected.
%
% see also: PKG.E_GETTFLIST, GLY.ENZONTO, SC_TFACTIVITY, GLY.ENRICH

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    opts.Species (1, 1) string = "hs"
    opts.GeneSets = []
    opts.NumBins (1, 1) double {mustBeGreaterThanOrEqual(opts.NumBins, 2)} = 5
    opts.MaxCells (1, 1) double {mustBePositive} = 3000
    opts.MinDetected (1, 1) double {mustBeInRange(opts.MinDetected, 0, 1)} = 0.02
    opts.Validate (1, 1) logical = true
    opts.Seed (1, 1) double = 42
end

if numel(genelist) ~= size(X, 1)
    error("GLY:TFMI:GeneCount", ...
        "GENELIST length (%d) must equal the number of rows of X (%d).", ...
        numel(genelist), size(X, 1));
end

TT = pkg.e_gettflist(opts.Species);
tfAll = unique(upper(string(TT.tf)));
glyAll = i_glycogenes(opts.GeneSets);

up = upper(genelist);
isTF = ismember(up, tfAll);
isGly = ismember(up, glyAll);
if ~any(isTF) || ~any(isGly)
    error("GLY:TFMI:NoOverlap", ...
        "Found %d TFs and %d glycogenes in GENELIST; need both.", ...
        sum(isTF), sum(isGly));
end

rng(opts.Seed);
nCells = size(X, 2);
cols = 1:nCells;
if nCells > opts.MaxCells
    cols = sort(randsample(nCells, opts.MaxCells))';
end

% Detection floor. An MI estimate over a gene seen in a handful of cells is
% noise with a number attached, and it would sit at the top of the ranking
% as readily as anywhere else.
Xs = full(double(X(:, cols)));
n = numel(cols);
frac = sum(Xs > 0, 2) / n;
keepTF = isTF & frac >= opts.MinDetected;
keepGly = isGly & frac >= opts.MinDetected;
if ~any(keepTF) || ~any(keepGly)
    error("GLY:TFMI:AllFiltered", ...
        ['Nothing survived MinDetected=%g: %d TFs and %d glycogenes ', ...
        'remain. Lower it or use a deeper dataset.'], ...
        opts.MinDetected, sum(keepTF), sum(keepGly));
end

tfs = genelist(keepTF);
glys = genelist(keepGly);
A = i_bin(Xs(keepTF, :)', opts.NumBins);      % cells x TFs
B = i_bin(Xs(keepGly, :)', opts.NumBins);     % cells x glycogenes

M = i_mi(A, B, opts.NumBins);

info = struct();
info.tfs = tfs;
info.glycogenes = glys;
info.nCells = n;
info.numBins = opts.NumBins;
info.fracDetected = struct('tf', frac(keepTF), 'glyco', frac(keepGly));

% Ranked summaries, so the matrix is usable without further work.
[~, ord] = sort(sum(M, 2), "descend");
info.topTFs = table(tfs(ord), sum(M(ord, :), 2), ...
    VariableNames = ["tf", "totalMI"]);
[v, li] = maxk(M(:), 25);
[ri, ci] = ind2sub(size(M), li);
info.topPairs = table(tfs(ri), glys(ci), v, ...
    VariableNames = ["tf", "glycogene", "mi"]);

if opts.Validate
    info.auprc = i_validate(M, Xs(keepTF, :), Xs(keepGly, :), tfs, glys, TT);
end

end


% ----------------------------------------------------------------------
function g = i_glycogenes(G)
if isempty(G)
    [~, ~, g] = gly.enzonto(Include = "all");
    g = upper(g);
    return;
end
if ~istable(G) || ~ismember("Genes", G.Properties.VariableNames)
    error("GLY:TFMI:BadGeneSets", ...
        "GeneSets must be a table with a Genes column.");
end
g = strings(0, 1);
for k = 1:height(G)
    v = strtrim(split(string(G.Genes(k)), ","));
    g = [g; v(strlength(v) > 0)]; %#ok<AGROW>
end
g = unique(upper(g));
end


function Bi = i_bin(V, nb)
% V is cells-by-genes. Bin 1 is "not detected"; the nonzero values go into
% nb-1 quantile bins. Equal-frequency binning over everything would put the
% cut points inside the zeros and waste the grid on them.
[n, p] = size(V);
Bi = ones(n, p);
for k = 1:p
    v = V(:, k);
    nz = v > 0;
    if ~any(nz), continue, end
    x = v(nz);
    edges = unique(quantile(x, linspace(0, 1, nb)));
    if numel(edges) < 2
        Bi(nz, k) = 2;
    else
        edges(1) = -Inf; edges(end) = Inf;
        b = discretize(x, edges);
        b(isnan(b)) = 1;
        Bi(nz, k) = 1 + b;
    end
end
Bi = min(Bi, nb);
end


function M = i_mi(A, B, nb)
% MI for every column pair, as a sum over the bin grid. For each bin pair
% (i,j) the joint count matrix is an indicator product, so the whole grid is
% nb^2 matrix multiplications rather than one histogram per gene pair --
% which for 1200 TFs by 400 glycogenes would be half a million of them.
n = size(A, 1);
nT = size(A, 2); nG = size(B, 2);
Ai = cell(nb, 1); Bj = cell(nb, 1);
pA = zeros(nT, nb); pB = zeros(nG, nb);
for i = 1:nb
    Ai{i} = single(A == i);
    Bj{i} = single(B == i);
    pA(:, i) = sum(Ai{i}, 1)' / n;
    pB(:, i) = sum(Bj{i}, 1)' / n;
end

M = zeros(nT, nG);
for i = 1:nb
    if ~any(pA(:, i)), continue, end
    for j = 1:nb
        if ~any(pB(:, j)), continue, end
        Pij = double(Ai{i}' * Bj{j}) / n;          % nT x nG joint
        denom = pA(:, i) * pB(:, j)';
        ok = Pij > 0 & denom > 0;
        if ~any(ok, "all"), continue, end
        term = zeros(nT, nG);
        term(ok) = Pij(ok) .* log(Pij(ok) ./ denom(ok));
        M = M + term;
    end
end
M = max(M, 0);      % clip the -0 that floating point leaves behind
end


function T = i_validate(M, Xtf, Xgly, tfs, glys, TT)
% Does MI recover known TF-target edges better than a correlation does, on
% THIS data? The comparison the paper makes, against DoRothEA rather than
% TFLink because that is the network this toolbox ships.
edges = unique(upper(string(TT.tf)) + "|" + upper(string(TT.target)));
[ti, gi] = ndgrid(1:numel(tfs), 1:numel(glys));
key = upper(tfs(ti(:))) + "|" + upper(glys(gi(:)));
y = ismember(key, edges);

if ~any(y) || all(y)
    T = table("none", NaN, NaN, ...
        VariableNames = ["estimator", "auprc", "prevalence"]);
    return;
end

P = corr(Xtf', Xgly');
S = corr(Xtf', Xgly', Type = "Spearman");
scores = {M(:), abs(P(:)), abs(S(:))};
names = ["mutual information"; "|Pearson|"; "|Spearman|"];
a = NaN(3, 1);
for k = 1:3
    a(k) = i_auprc(scores{k}, y);
end
prev = mean(y);
T = table(names, a, repmat(prev, 3, 1), ...
    VariableNames = ["estimator", "auprc", "prevalence"]);
T = sortrows(T, "auprc", "descend");
end


function a = i_auprc(score, y)
% Average precision: the area under the precision-recall curve as the mean
% precision at each true positive. Ties are broken arbitrarily by SORT,
% which is the usual convention and immaterial at this scale.
ok = isfinite(score);
score = score(ok); y = y(ok);
[~, ord] = sort(score, "descend");
y = y(ord);
tp = cumsum(y);
prec = tp ./ (1:numel(y))';
a = sum(prec(y)) / max(sum(y), 1);
end
