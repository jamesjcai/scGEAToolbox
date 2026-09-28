function [T, info] = ml_GSEA(X, genelist, grp, setmatrx, setnames, setgenes, opts)
%ML_GSEA  Standard GSEA with sample-label (phenotype) permutation.
%
%   T = RUN.ML_GSEA(X, genelist, grp, setmatrx, setnames, setgenes) runs the
%   original Broad GSEA on a genes-by-SAMPLES matrix: rank every gene by how
%   well it separates the two phenotypes, score each gene set with the
%   weighted Kolmogorov-Smirnov running sum, and build the null by permuting
%   the sample labels and repeating the whole thing.
%
%   WHY THE LABEL PERMUTATION MATTERS. SC_GSETTEST and SC_FGSEA permute gene
%   labels, which assumes genes are independent. They are not: a gene set is
%   a set precisely because its genes are co-expressed, so the null is too
%   narrow and the false discovery rate is understated. Permuting sample
%   labels instead leaves every gene-gene correlation exactly where it was
%   and asks the question the experiment actually poses - would this
%   separation appear if the phenotype labels were meaningless? That is the
%   right null whenever the data support it.
%
%   WHAT IT COSTS. It needs replicate samples, because the null is built
%   from relabellings of them: with n per group there are only
%   nchoosek(2n, n) distinct assignments, and below about 7 per group the
%   null is too coarse to resolve a small p-value at all. The function warns
%   when a group is smaller than that, and again when NumPerm exceeds the
%   number of distinct relabellings.
%
%   SINGLE-CELL DATA IS NOT SAMPLE-LEVEL DATA. Cells from one donor are not
%   independent replicates of that donor's phenotype. Passing cells as
%   columns treats thousands of correlated measurements as thousands of
%   samples and will call nearly every set significant. Aggregate to
%   pseudobulk first - one column per donor, sample or well - and pass that.
%   The function warns above 200 columns for this reason.
%
%   USAGE:
%     % pseudobulk: one column per donor, two phenotypes
%     [setmatrx, setnames, setgenes] = pkg.e_getgenesets('MSIGDB', 'human');
%     T = run.ml_GSEA(P, genelist, isTreated, setmatrx, setnames, setgenes);
%
%     % an ad-hoc collection, passed as gene lists
%     T = run.ml_GSEA(P, g, grp, {["A" "B"], ["C" "D"]}, ["s1" "s2"], []);
%
%   INPUTS:
%     X        - genes-by-samples expression on a log scale (log-CPM, VST,
%                log1p of normalised counts). Not raw counts: the ranking
%                metric is a difference of means over a spread, and on the
%                raw scale that is dominated by the most highly expressed
%                genes. May be sparse.
%     genelist - G-by-1 gene symbols, one per row of X.
%     grp      - S-by-1 phenotype, one entry per column of X. Either a
%                logical vector (true = the class of interest) or any
%                two-level grouping (string, cellstr, categorical, numeric).
%     setmatrx - nSets-by-nSetGenes membership matrix from
%                PKG.E_GETGENESETS, or an nSets-by-1 cell array of gene-name
%                lists (SETGENES is then ignored).
%     setnames - nSets-by-1 set names ([] to auto-generate).
%     setgenes - nSetGenes-by-1 gene symbols for the columns of SETMATRX.
%
%   NAME-VALUE ARGUMENTS:
%     PositiveClass - which level of GRP is the class of interest, i.e. the
%                     one a positive ES points to. Default is the second of
%                     the two levels in sorted order, which puts "treated"
%                     above "control" and "case" above "baseline". Ignored
%                     when GRP is logical, where true is the class.
%     Metric        - gene ranking statistic:
%                     "s2n" (default) signal-to-noise, (m1-m2)/(s1+s2) with
%                       the Broad's standard-deviation floor, which stops a
%                       gene with almost no spread from ranking first on a
%                       trivial difference. This is GSEA's own default.
%                     "ttest" Welch t-statistic; sharper when the group
%                       variances genuinely differ.
%                     "diff" plain difference of means, with no variance
%                       stabilisation at all.
%     Weight        - exponent p on |metric| in the running sum (default 1).
%                     0 gives the classic unweighted KS statistic.
%     MinSize       - minimum genes of a set present in GENELIST (default
%                     15, the GSEA convention).
%     MaxSize       - maximum (default 500). Very large sets otherwise crowd
%                     the top of the result and say little.
%     NumPerm       - sample-label permutations (default 1000). Unlike the
%                     gene-label nulls of SC_GSETTEST this cannot be raised
%                     without limit: the distinct relabellings run out.
%     Sort          - sort rows by ascending PValue (default true).
%     Verbose       - print a short summary (default true).
%
%   OUTPUTS:
%     T    - one row per tested set:
%              SetName     set name
%              SetSize     genes of the set found in GENELIST
%              ES          enrichment score; positive means the set sits
%                          toward the PositiveClass end of the ranking
%              NES         ES divided by the mean of the same-sign null,
%                          which is what makes sets of different size
%                          comparable
%              PValue      permutation p-value against the same-sign null
%              FDR         Benjamini-Hochberg adjusted p-value
%              LeadingEdge up to ten set genes at the ES peak, comma
%                          separated; INFO holds all of them
%     info - struct with the membership matrix, the gene universe, the
%            observed metric and ranking, the full leading edges, and the
%            options used. Its LeadingEdge follows the rows of T.
%
%   THE P-VALUE FLOOR. Only the half of the null lying on the observed
%   side counts, so with NumPerm draws the smallest reportable p-value is
%   near 2/NumPerm, not 1/NumPerm - and unlike the gene-permutation path
%   there is no tail extrapolation here, because fitting a tail needs
%   thousands of independent draws and the sample labels cannot supply
%   them. Sets at the floor are tied, and the ranking among them is not
%   resolved. If that resolution is what is wanted, and the independence
%   assumption is acceptable, use SC_FGSEA or
%   SC_GSETTEST(Method="gsea", PValueTail="gpd") instead.
%
% REF: Subramanian et al. (2005) PNAS 102:15545; Mootha et al. (2003)
%      Nat Genet 34:267; Zyla et al. (2017) Bioinformatics 33:5382 on the
%      choice of ranking metric.
%
% See also SC_GSETTEST, SC_FGSEA, PKG.E_GETGENESETS, PKG.E_GSEASCORE.

arguments
    X {mustBeNumeric, mustBeNonempty}
    genelist (:, 1) string
    grp
    setmatrx
    setnames = []
    setgenes = []
    opts.PositiveClass = []
    opts.Metric (1, 1) string {mustBeMember(opts.Metric, ...
        ["s2n", "ttest", "diff"])} = "s2n"
    opts.Weight (1, 1) double {mustBeNonnegative} = 1
    opts.MinSize (1, 1) double {mustBePositive} = 15
    opts.MaxSize (1, 1) double {mustBePositive} = 500
    opts.NumPerm (1, 1) double {mustBePositive, mustBeInteger} = 1000
    opts.Sort (1, 1) logical = true
    opts.Verbose (1, 1) logical = true
end

if size(X, 1) ~= numel(genelist)
    error("run:ml_GSEA:SizeMismatch", ...
        "X has %d rows and GENELIST has %d entries. X must be " + ...
        "genes-by-samples.", size(X, 1), numel(genelist));
end

numSample = size(X, 2);
[isPos, classNames] = i_resolvephenotype(grp, numSample, opts.PositiveClass);
n1 = nnz(isPos);
n2 = numSample - n1;

if min(n1, n2) < 2
    error("run:ml_GSEA:TooFewSamples", ...
        "Group ""%s"" has %d sample(s) and ""%s"" has %d. Label " + ...
        "permutation needs at least 2 in each, and realistically 7.", ...
        classNames(2), n1, classNames(1), n2);
end

% Through GAMMALN rather than NCHOOSEK, which warns about precision once
% the count passes flintmax - and it does, for any decently sized design.
% Only its order of magnitude is ever used here.
numDistinct = round(exp(gammaln(numSample + 1) - gammaln(n1 + 1) - ...
    gammaln(n2 + 1)));
if min(n1, n2) < 7
    warning("run:ml_GSEA:FewReplicates", ...
        "Only %d sample(s) in the smaller group. The %g distinct " + ...
        "relabellings cannot resolve a p-value below %.3g, so little " + ...
        "here will survive an FDR cutoff on a large collection. " + ...
        "SC_FGSEA on a differential statistic is the alternative.", ...
        min(n1, n2), numDistinct, 1/(numDistinct + 1));
end
if opts.NumPerm > numDistinct
    % Drawing more labellings than exist does not add information; it just
    % resamples the same ones, and the p-value stops improving.
    warning("run:ml_GSEA:PermutationsExhausted", ...
        "NumPerm is %d but only %g distinct relabellings exist. The " + ...
        "extra draws are repeats and buy no resolution.", ...
        opts.NumPerm, numDistinct);
end
if numSample > 200
    % The columns are meant to be biological replicates. Hundreds of them
    % almost always means cells were passed instead of pseudobulk, and the
    % test then treats correlated cells as independent samples.
    warning("run:ml_GSEA:ManySamples", ...
        "X has %d columns. If those are cells rather than samples, " + ...
        "aggregate to pseudobulk first: cells from one donor are not " + ...
        "replicates of that donor, and every set will look significant.", ...
        numSample);
end

% ---- gene universe --------------------------------------------------
genelist = upper(genelist);
keep = strlength(genelist) > 0;
X = X(keep, :);
genelist = genelist(keep);

[ug, ~, ic] = unique(genelist);
if numel(ug) < numel(genelist)
    % The ranking is recomputed once per permutation, so the representative
    % copy has to be chosen once and stay chosen - "most extreme statistic",
    % the way SC_GSETTEST picks it, would change from draw to draw. Mean
    % expression is a property of the gene rather than of the labelling.
    warning("run:ml_GSEA:DuplicateGenes", ...
        "%d duplicate gene symbols collapsed to the copy with the " + ...
        "highest mean expression.", numel(genelist) - numel(ug));
    rowMean = full(mean(X, 2));
    pick = accumarray(ic, (1:numel(ic))', [], @(v) i_argmax(rowMean, v));
    X = X(pick, :);
    genelist = ug;
end

numGene = numel(genelist);
if numGene < 10
    error("run:ml_GSEA:TooFewGenes", ...
        "Only %d usable genes; an enrichment score is meaningless.", numGene);
end

% ---- membership over that universe ----------------------------------
[setmatrx, setnames, setgenes] = pkg.i_normalizegenesets(setmatrx, ...
    setnames, setgenes);
if ~issparse(setmatrx)
    setmatrx = sparse(double(setmatrx));
end
[tf, loc] = ismember(genelist, upper(setgenes));
M = sparse(size(setmatrx, 1), numGene);
M(:, tf) = setmatrx(:, loc(tf));

msize = full(sum(M ~= 0, 2));
ok = msize >= opts.MinSize & msize <= min(opts.MaxSize, numGene - 1);
if ~any(ok)
    error("run:ml_GSEA:NoSets", ...
        "No gene set has between %g and %g of its genes in GENELIST " + ...
        "(largest overlap was %d). Check that GENELIST and SETGENES use " + ...
        "the same symbol namespace.", opts.MinSize, opts.MaxSize, max(msize));
end
M = spones(M(ok, :));
names = setnames(ok);
m = full(sum(M, 2));
numSet = numel(m);
offset = [0; cumsum(m)];

% ---- per-sample sufficient statistics --------------------------------
% Every permutation needs a group mean and spread, and both are sums over
% the columns of one group. Holding X and X.^2 and multiplying by the group
% indicator turns each relabelling into two matrix-vector products, so
% nothing is re-extracted or densified per draw.
X = double(X);
Xsq = X.^2;
rowSum = full(sum(X, 2));
rowSumSq = full(sum(Xsq, 2));

[stat, ord, w] = i_rank(X, Xsq, rowSum, rowSumSq, double(isPos), n1, n2, opts);

% ---- observed scores and leading edges -------------------------------
qObs = i_positions(M, ord, numGene);
es = zeros(numSet, 1);
peak = zeros(numSet, 1);
for k = 1:numSet
    [es(k), peak(k)] = pkg.e_gseascore(w, ...
        qObs(offset(k)+1:offset(k+1)), numGene, m(k));
end

rankedGenes = genelist(ord);
leadingEdge = cell(numSet, 1);
for k = 1:numSet
    q = qObs(offset(k)+1:offset(k+1));
    if es(k) > 0
        leadingEdge{k} = rankedGenes(q(q <= peak(k)));
    elseif es(k) < 0
        leadingEdge{k} = rankedGenes(q(q >= peak(k)));
    else
        leadingEdge{k} = strings(0, 1);
    end
end

% ---- null by relabelling ---------------------------------------------
% Accumulated rather than stored: the same-sign mean and the exceedance
% count are all the null is asked for, and a NumPerm-by-nSets array over a
% large collection would be hundreds of megabytes of nothing else.
side = sign(es);
side(side == 0) = 1;
posSum = zeros(numSet, 1);
posCount = zeros(numSet, 1);
negSum = zeros(numSet, 1);
negCount = zeros(numSet, 1);
numHit = zeros(numSet, 1);
esPerm = zeros(numSet, 1);

if opts.Verbose
    fprintf("ml_GSEA: %d sets over %d samples, %d label permutations\n", ...
        numSet, numSample, opts.NumPerm);
end

for b = 1:opts.NumPerm
    gperm = double(isPos(randperm(numSample)));
    [~, ordp, wp] = i_rank(X, Xsq, rowSum, rowSumSq, gperm, n1, n2, opts);
    qp = i_positions(M, ordp, numGene);
    for k = 1:numSet
        esPerm(k) = pkg.e_gseascore(wp, ...
            qp(offset(k)+1:offset(k+1)), numGene, m(k));
    end
    isUpPerm = esPerm > 0;
    isDownPerm = esPerm < 0;
    posSum = posSum + esPerm.*isUpPerm;
    posCount = posCount + isUpPerm;
    negSum = negSum + esPerm.*isDownPerm;
    negCount = negCount + isDownPerm;
    numHit = numHit + (side.*esPerm >= side.*es);
end

% ---- NES, p-value, FDR ------------------------------------------------
% Both are read off the half of the null lying on the observed side, which
% is how GSEA normalises: the positive and negative halves have different
% magnitudes, and pooling them would scale every score by a mixture that
% depends on how the null happened to split.
isUp = es >= 0;
refMean = nan(numSet, 1);
refCount = posCount;
refCount(~isUp) = negCount(~isUp);
refMean(isUp & posCount > 0) = posSum(isUp & posCount > 0) ./ ...
    posCount(isUp & posCount > 0);
refMean(~isUp & negCount > 0) = negSum(~isUp & negCount > 0) ./ ...
    negCount(~isUp & negCount > 0);

nes = nan(numSet, 1);
usable = isfinite(refMean) & refMean ~= 0;
nes(usable) = es(usable) ./ abs(refMean(usable));

p = (1 + numHit) ./ (1 + refCount);
fdr = pkg.e_fdr(p);

leadingEdgeText = strings(numSet, 1);
for k = 1:numSet
    g = leadingEdge{k};
    if ~isempty(g)
        leadingEdgeText(k) = strjoin(g(1:min(10, numel(g))), ",");
    end
end

T = table(string(names(:)), m, es, nes, p, fdr(:), leadingEdgeText, ...
    VariableNames=["SetName", "SetSize", "ES", "NES", "PValue", "FDR", ...
    "LeadingEdge"]);
if opts.Sort
    [T, sortIdx] = sortrows(T, "PValue");
    leadingEdge = leadingEdge(sortIdx);
end

if opts.Verbose
    fprintf("  metric=%s, ""%s"" (%d) vs ""%s"" (%d)\n", opts.Metric, ...
        classNames(2), n1, classNames(1), n2);
    fprintf("  universe %d genes; %d of %d sets tested (size %d-%d)\n", ...
        numGene, numSet, numel(ok), min(m), max(m));
    fprintf("  %d sets at FDR<0.05, %d at FDR<0.25\n", ...
        sum(T.FDR < 0.05), sum(T.FDR < 0.25));
end

if nargout > 1
    info = struct(M=M, Universe=genelist, Stat=stat, Ranking=rankedGenes, ...
        SetIndex=find(ok), LeadingEdge={leadingEdge}, Phenotype=isPos, ...
        Classes=classNames, Options=opts);
end

end

% =========================================================================
function idx = i_argmax(v, group)
% The member of GROUP holding its largest value in V.
[~, k] = max(v(group));
idx = group(k);
end

% =========================================================================
function [isPos, classNames] = i_resolvephenotype(grp, numSample, positiveClass)
% Reduce any two-level grouping to a logical, and record which level is
% which so that the sign of ES can be reported in the user's own terms.
if numel(grp) ~= numSample
    error("run:ml_GSEA:SizeMismatch", ...
        "GRP has %d entries but X has %d columns.", numel(grp), numSample);
end

if islogical(grp)
    isPos = grp(:);
    classNames = ["false", "true"];
else
    g = string(grp(:));
    levels = unique(g);
    levels = levels(strlength(levels) > 0 & ~ismissing(levels));
    if numel(levels) ~= 2
        error("run:ml_GSEA:NotTwoGroups", ...
            "GRP has %d distinct levels (%s). A phenotype permutation " + ...
            "test compares exactly two.", numel(levels), ...
            strjoin(levels(1:min(5, numel(levels))), ", "));
    end
    if isempty(positiveClass)
        positive = levels(2);
    else
        positive = string(positiveClass);
        if ~ismember(positive, levels)
            error("run:ml_GSEA:BadPositiveClass", ...
                "PositiveClass ""%s"" is not a level of GRP (%s).", ...
                positive, strjoin(levels, ", "));
        end
    end
    isPos = g == positive;
    classNames = [levels(levels ~= positive), positive];
end

if all(isPos) || ~any(isPos)
    error("run:ml_GSEA:NotTwoGroups", ...
        "Every sample is in the same group; there is nothing to compare.");
end
end

% =========================================================================
function [stat, ord, w] = i_rank(X, Xsq, rowSum, rowSumSq, g, n1, n2, opts)
% Rank the genes by how well they separate the two groups under labelling G.
sum1 = full(X*g);
sumSq1 = full(Xsq*g);
sum2 = rowSum - sum1;
sumSq2 = rowSumSq - sumSq1;
mean1 = sum1/n1;
mean2 = sum2/n2;

switch opts.Metric
    case "diff"
        stat = mean1 - mean2;

    case "ttest"
        var1 = max((sumSq1 - n1*mean1.^2)/(n1 - 1), 0);
        var2 = max((sumSq2 - n2*mean2.^2)/(n2 - 1), 0);
        se = sqrt(var1/n1 + var2/n2);
        stat = (mean1 - mean2) ./ se;
        stat(se == 0) = 0;          % constant in both groups: no signal

    case "s2n"
        sd1 = sqrt(max((sumSq1 - n1*mean1.^2)/(n1 - 1), 0));
        sd2 = sqrt(max((sumSq2 - n2*mean2.^2)/(n2 - 1), 0));
        % The Broad's standard-deviation floor. Without it a gene that
        % happens to be nearly constant in both groups divides a tiny
        % difference by a near-zero spread and tops the ranking on noise.
        sd1 = max(sd1, 0.2*abs(mean1));
        sd1(sd1 == 0) = 0.2;
        sd2 = max(sd2, 0.2*abs(mean2));
        sd2(sd2 == 0) = 0.2;
        stat = (mean1 - mean2) ./ (sd1 + sd2);

    otherwise
        % Unreachable: MUSTBEMEMBER in the arguments block admits no other
        % value. Kept so that adding a metric there without adding it here
        % fails loudly rather than ranking by whatever was left over.
        error("run:ml_GSEA:UnknownMetric", ...
            "Metric ""%s"" has no implementation.", opts.Metric);
end

[ssort, ord] = sort(stat, "descend");
w = abs(ssort).^opts.Weight;
end

% =========================================================================
function q = i_positions(M, ord, numGene)
% Ranks held by each set's genes, ascending within a set and concatenated
% set by set. Transposing the reordered membership matrix puts each set in
% one column, and FIND walks a sparse matrix column by column with the row
% indices already in order - so this replaces one sort per set per
% permutation with a single pass over the whole collection.
idx = find(M(:, ord).');
q = rem(idx - 1, numGene) + 1;
end
