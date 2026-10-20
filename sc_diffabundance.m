function [T, Counts] = sc_diffabundance(celltype, group, sampleid, methodid, ...
    nperm, transform)
% SC_DIFFABUNDANCE  Differential cell-type abundance testing between conditions.
%
% Detects cell types whose relative abundance shifts between two (or more)
% conditions/regions - e.g., a lineage that is specifically depleted in one
% region while other cell types are unchanged, the hallmark of cell
% competition (Jiang et al., gliomagenesis, bioRxiv 2026.01.15.699808).
%
% USAGE:
%   T = sc_diffabundance(celltype, group);
%   T = sc_diffabundance(celltype, group, sampleid);         % replicate-aware
%   [T, Counts] = sc_diffabundance(celltype, group, sampleid, methodid, nperm);
%
% INPUTS:
%   celltype  - C-by-1 per-cell cell-type labels (string, cellstr, or categorical)
%   group     - C-by-1 per-cell condition/region labels (e.g., "core"/"infiltrating").
%               The first group (sorted) is the reference for fold-change.
%   sampleid  - (optional) C-by-1 replicate/sample labels. When provided with
%               methodid 1, 4 or 5, testing is done across samples (recommended),
%               which properly accounts for biological replication.
%   methodid  - test to use (default 1):
%               1 = auto: propeller (5) if sampleid has >=2 samples/group,
%                   otherwise cell-level permutation (2)
%               2 = cell-level label-permutation test (no replicates)
%               3 = Fisher's exact test per cell type (two groups only)
%               4 = replicate test on per-unit proportions, where a unit
%                   is one (sample, group) combination. For two groups:
%                   Wilcoxon signed-rank when every sample appears in both
%                   groups (a paired design, e.g. one tumour dissected into
%                   core and infiltrating), Wilcoxon rank-sum otherwise.
%                   Kruskal-Wallis for >2 groups. Requires sampleid.
%               5 = propeller (Phipson et al. 2022): per-unit proportions
%                   are transformed (logit or arcsin square root), fitted
%                   by a linear model per cell type, and tested with an
%                   empirical-Bayes moderated t (two groups) or moderated
%                   F (>2 groups) that borrows variance across cell types.
%                   Units are (sample, group) combinations as in 4; when
%                   every sample is in every group the model also removes
%                   sample effects (paired / block design). Unlike the
%                   rank tests in 4 it has power with 2-3 replicates per
%                   group. Matches speckle::propeller(robust = FALSE);
%                   speckle's default robust = TRUE is not implemented.
%                   Requires sampleid.
%   nperm     - number of permutations for methodid 2 (default 2000)
%   transform - methodid 5 only: "logit" (default; speckle's default) or
%               "asin". The paper finds asin less sensitive to outlier
%               samples.
%
% OUTPUTS:
%   T       - table with one row per cell type:
%               CellType, Prop_<g1>..Prop_<gk> (mean proportions per group),
%               Log2FC (two-group only; log2 of prop ratio vs reference),
%               Pvalue, FDR, Direction ("enriched"/"depleted"/"n.s.")
%             Rows are sorted by FDR ascending.
%   Counts  - cell-type-by-group table of raw cell counts (cross-tabulation).
%
%             T.Properties.Description names the test that was performed.
%
% see also: SC_COMPETITIONSCORE, SC_CELLSCORE, PKG.E_FDR_BH

if nargin < 2, error("USAGE: T = sc_diffabundance(celltype, group);"); end
if nargin < 3, sampleid = []; end
if nargin < 4 || isempty(methodid), methodid = 1; end
if nargin < 5 || isempty(nperm), nperm = 2000; end
if nargin < 6 || isempty(transform), transform = "logit"; end
transform = lower(string(transform));

celltype = string(celltype(:));
group = string(group(:));
if numel(celltype) ~= numel(group)
    error("celltype and group must have the same number of cells.");
end
if ~isempty(sampleid)
    sampleid = string(sampleid(:));
    if numel(sampleid) ~= numel(celltype)
        error("sampleid must have the same number of cells as celltype.");
    end
end

ctypes = unique(celltype);
grps = unique(group);
nCt = numel(ctypes);
nG = numel(grps);
if nG < 2
    error("group must contain at least two distinct conditions.");
end

% Raw cross-tabulation (cell types x groups) for reporting.
countMat = zeros(nCt, nG);
for j = 1:nG
    inGrp = group == grps(j);
    countMat(:, j) = countByType(celltype(inGrp), ctypes);
end
Counts = array2table(countMat, RowNames=cellstr(ctypes), ...
    VariableNames=matlab.lang.makeValidName("n_"+grps));

% Overall proportion of each cell type within each group.
propMat = countMat ./ sum(countMat, 1);

% Resolve "auto" method.
if methodid == 1
    if ~isempty(sampleid) && hasReplicates(sampleid, group, grps)
        methodid = 5;
    else
        methodid = 2;
    end
end

switch methodid
    case 2
        pvals = i_permtest(celltype, group, ctypes, grps, nperm);
        testname = "cell-level label permutation";
    case 3
        if nG ~= 2
            error("Fisher's exact test (methodid 3) requires exactly two groups.");
        end
        pvals = i_fisher(countMat);
        testname = "Fisher's exact";
    case 4
        if isempty(sampleid)
            error("methodid 4 requires sampleid for replicate-aware testing.");
        end
        [pvals, testname] = i_replicatetest(celltype, group, sampleid, ...
            ctypes, grps);
    case 5
        if isempty(sampleid)
            error("methodid 5 (propeller) requires sampleid.");
        end
        [pvals, testname, effect] = i_propeller(celltype, group, ...
            sampleid, ctypes, grps, transform);
    otherwise
        error("Unknown methodid %d. Use 1, 2, 3, 4, or 5.", methodid);
end

fdr = pkg.e_fdr(pvals);
fdr = fdr(:);

% Assemble output table.
T = table(ctypes, VariableNames="CellType");
for j = 1:nG
    T.("Prop_"+matlab.lang.makeValidName(grps(j))) = propMat(:, j);
end
if nG == 2
    log2fc = log2((propMat(:, 2) + eps) ./ (propMat(:, 1) + eps));
    T.Log2FC = log2fc;
else
    log2fc = nan(nCt, 1);
end
T.Pvalue = pvals(:);
T.FDR = fdr;

direction = repmat("n.s.", nCt, 1);
sig = fdr < 0.05;
if nG == 2
    % propeller tests sample-level proportions, which a pooled Log2FC
    % dominated by one large sample can contradict; take its own sign.
    if methodid == 5
        shift = effect;
    else
        shift = log2fc;
    end
    direction(sig & shift > 0) = "enriched";
    direction(sig & shift < 0) = "depleted";
else
    direction(sig) = "changed";
end
T.Direction = direction;

[~, ord] = sort(T.FDR, "ascend");
T = T(ord, :);
T.Properties.Description = sprintf("sc_diffabundance methodid %d: %s", ...
    methodid, testname);

end


%% ---- local helpers ----
function counts = countByType(labels, ctypes)
% Count occurrences of each cell type (in ctypes order).
[tf, loc] = ismember(labels, ctypes);
counts = accumarray(loc(tf), 1, [numel(ctypes), 1]);
end


function tf = hasReplicates(sampleid, group, grps)
% True when every group contains at least two distinct samples.
tf = true;
for j = 1:numel(grps)
    if numel(unique(sampleid(group == grps(j)))) < 2
        tf = false;
        return;
    end
end
end


function pvals = i_permtest(celltype, group, ctypes, grps, nperm)
% Cell-level label-permutation test. Statistic per cell type is the
% two-group proportion difference, or a chi-square-like deviation for >2
% groups. The null is built by shuffling group labels across all cells.
nCells = numel(celltype);
nCt = numel(ctypes);
nG = numel(grps);

[~, ctIdx] = ismember(celltype, ctypes);
O = sparse(1:nCells, ctIdx, 1, nCells, nCt); % cells x celltypes one-hot

[~, gIdx] = ismember(group, grps);
obsStat = i_stat(O, gIdx, nG);

nullGE = zeros(nCt, 1); % count of null stats >= observed
for p = 1:nperm
    permG = gIdx(randperm(nCells));
    nullStat = i_stat(O, permG, nG);
    nullGE = nullGE + (nullStat >= obsStat - 1e-12);
end
pvals = (1 + nullGE) / (1 + nperm);
end


function stat = i_stat(O, gIdx, nG)
% Per-cell-type test statistic given a one-hot cell-type matrix O
% (cells x celltypes) and a group-index vector gIdx.
nCells = size(O, 1);
G = sparse(1:nCells, gIdx, 1, nCells, nG); % cells x groups one-hot
counts = full(O' * G); % celltypes x groups
groupTot = sum(counts, 1); % 1 x groups
props = counts ./ groupTot; % proportion of each celltype within each group
if nG == 2
    stat = abs(props(:, 2) - props(:, 1));
else
    expected = mean(props, 2);
    stat = sum((props - expected).^2 ./ (expected + eps), 2);
end
end


function pvals = i_fisher(countMat)
% Fisher's exact test per cell type on a 2x2 (in-type vs rest) x group table.
nCt = size(countMat, 1);
grpTot = sum(countMat, 1);
pvals = ones(nCt, 1);
for k = 1:nCt
    inType = countMat(k, :);
    rest = grpTot - inType;
    tbl = [inType(1), inType(2); rest(1), rest(2)];
    [~, pvals(k)] = fishertest(tbl);
end
end


function [pvals, testname] = i_replicatetest(celltype, group, sampleid, ...
    ctypes, grps)
% Replicate-aware test on per-unit cell-type proportions.
%
%   A replicate unit is one (SAMPLE, GROUP) combination, not one sample.
%
%   This used to assume the opposite:
%
%       gi = unique(group(inS));
%       sampleGroup(i) = gi(1); % a sample belongs to a single group
%
%   taking the first of however many groups the sample's cells span, with a
%   comment stating the assumption instead of checking it. In the design
%   this function's own header cites -- one tumour dissected into core and
%   infiltrating regions -- every sample spans both groups, so every sample
%   was labelled with the alphabetically first group, the second group's
%   set came back empty, and the isempty guard below turned that into
%   pvals(k) = 1 for every cell type. The caller got a full table reading
%   Pvalue 1, FDR 1 and Direction "n.s." throughout, with no warning that
%   no test had been run. Measured on three gliomas with one lineage 16x
%   depleted in the core (Log2FC 4.01): method 4 reported n.s. for all four
%   lineages while method 2 on the same cells found every one of them at
%   FDR < 0.01.
%
%   Splitting by (sample, group) also fixes a second problem the same line
%   caused: the proportions themselves pooled each sample's cells across
%   both groups, so the numbers fed to the test were region-blind even
%   before the labelling went wrong.
%
%   hasReplicates, which routes calls here, is satisfied -- indeed
%   maximally satisfied -- by samples appearing in both groups, so this was
%   the default path for exactly the design it could not handle.

nCt = numel(ctypes);
nG = numel(grps);

[countByUnit, ugIdx, unitSample, isBlocked] = i_units(celltype, group, ...
    sampleid, ctypes, grps);
nU = size(countByUnit, 1);
propByUnit = countByUnit ./ max(sum(countByUnit, 2), 1);
nPairs = numel(unique(unitSample));
isPaired = nG == 2 && isBlocked;

pvals = ones(nCt, 1);
if isPaired
    % Every sample appears in both groups: pair the units by sample.
    testname = "Wilcoxon signed-rank (paired by sample)";
    [~, ordA] = sort(unitSample(ugIdx == 1));
    [~, ordB] = sort(unitSample(ugIdx == 2));
    propA = propByUnit(ugIdx == 1, :);
    propB = propByUnit(ugIdx == 2, :);
    propA = propA(ordA, :);
    propB = propB(ordB, :);
    for k = 1:nCt
        if any(propB(:, k) ~= propA(:, k))
            pvals(k) = signrank(propA(:, k), propB(:, k));
        end
        % else: identical in every pair, so 1 is the right answer.
    end
    minP = 2^(1 - nPairs);
elseif nG == 2
    testname = "Wilcoxon rank-sum";
    for k = 1:nCt
        pvals(k) = ranksum(propByUnit(ugIdx == 1, k), ...
            propByUnit(ugIdx == 2, k));
    end
    n1 = nnz(ugIdx == 1);
    n2 = nnz(ugIdx == 2);
    minP = 2 / nchoosek(n1 + n2, n1);
else
    testname = "Kruskal-Wallis";
    for k = 1:nCt
        pvals(k) = kruskalwallis(propByUnit(:, k), ugIdx, "off");
    end
    minP = 0;   % not worth deriving for the >2-group case
end

% A rank test on few replicates has a floor on the p-value it can return,
% whatever the effect size. Saying so is the difference between "no
% difference" and "too few replicates to tell" -- and this function
% returning n.s. everywhere is exactly what the defect above looked like.
if minP > 0.05
    warning("sc_diffabundance:underpowered", ...
        ['%d replicate unit(s) with a %s cannot produce a two-sided ', ...
        'p-value below %.3g, so no cell type can reach FDR < 0.05 ', ...
        'however large its shift. Add samples, or use methodid 2 for a ', ...
        'cell-level test.'], nU, testname, minP);
end
end


function [countByUnit, ugIdx, unitSample, isBlocked] = i_units(celltype, ...
    group, sampleid, ctypes, grps)
% Cell-type counts per replicate unit, one unit per (sample, group)
% combination that actually contains cells. ISBLOCKED is true when there
% are at least two samples and every sample has exactly one unit in every
% group -- a paired (two groups) or randomized-block (more) design.
nCt = numel(ctypes);
nG = numel(grps);

% CHAR(31) is the ASCII unit separator, so it cannot occur in a label.
unitKey = sampleid + string(char(31)) + group;
[~, firstOf, unitOf] = unique(unitKey);
nU = numel(firstOf);

countByUnit = zeros(nU, nCt);
for i = 1:nU
    countByUnit(i, :) = countByType(celltype(unitOf == i), ctypes)';
end
unitGroup = group(firstOf);
unitSample = sampleid(firstOf);
[~, ugIdx] = ismember(unitGroup, grps);

% GRPS comes from GROUP, so every group has at least one unit.
if any(~ismember(1:nG, ugIdx))
    error("sc_diffabundance:groupWithoutUnits", ...
        "Internal error: a group ended up with no replicate unit.");
end

samples = unique(unitSample);
nS = numel(samples);
isBlocked = nS >= 2 && nU == nG*nS && ...
    all(arrayfun(@(s) isequal(sort(unitGroup(unitSample == s)), grps(:)), ...
    samples));
end


function [pvals, testname, effect] = i_propeller(celltype, group, sampleid, ...
    ctypes, grps, transform)
% propeller (Phipson et al. 2022, Bioinformatics 38:4720): transform the
% per-unit cell-type proportions, fit a linear model per cell type, and
% test the group effect with an empirical-Bayes moderated t (two groups)
% or moderated F (more groups) that borrows variance across cell types
% (Smyth 2004). Written from the two papers; speckle/limma were used only
% as a black box to check the output (tests/diffAbundanceTest.m).
%
% EFFECT is the estimated difference of the last group from the first on
% the transformed scale (two groups only), so its sign gives the direction.
nG = numel(grps);
[countByUnit, ugIdx, unitSample, isBlocked] = i_units(celltype, group, ...
    sampleid, ctypes, grps);
nU = size(countByUnit, 1);

switch transform
    case "logit"
        % Offset of 0.5 on the counts so that a zero count has a logit.
        p = (countByUnit + 0.5) ./ sum(countByUnit + 0.5, 2);
        Z = log(p ./ (1 - p));
    case "asin"
        p = countByUnit ./ sum(countByUnit, 2);
        Z = asin(sqrt(p));
    otherwise
        error("sc_diffabundance:badTransform", ...
            "transform must be ""logit"" or ""asin"", not ""%s"".", transform);
end

% Design: one mean per group. In a paired or block design every sample
% is in every group, so sample effects are estimable and are removed.
X = double(ugIdx(:) == 1:nG);
testname = sprintf("propeller (%s, moderated ", transform);
if isBlocked
    [~, ~, sIdx] = unique(unitSample);
    S = double(sIdx(:) == 2:max(sIdx));
    X = [X, S];
    blockNote = ", blocked by sample)";
else
    blockNote = ")";
end
nCoef = size(X, 2);
dfResid = nU - rank(X);
if dfResid < 1
    error("sc_diffabundance:noResidualDf", ...
        ['methodid 5 needs more replicate units than model coefficients ', ...
        '(%d units, %d coefficients). Add samples.'], nU, nCoef);
end

B = X \ Z;
resid = Z - X*B;
s2 = sum(resid.^2, 1)' / dfResid;

% Contrasts: every group against the first. The moderated F below does
% not depend on which basis of group differences is used.
C = [-ones(1, nG - 1); eye(nG - 1); zeros(nCoef - nG, nG - 1)];
est = (C'*B)';                               % celltypes x contrasts
covC = C'*((X'*X)\C);
unscaled = sqrt(diag(covC))';                % 1 x contrasts

[s2post, dfTotal, dfPost] = i_squeezevar(s2, dfResid);
tmod = est ./ (sqrt(s2post) .* unscaled);

if nG == 2
    testname = testname + "t" + blockNote;
    pvals = 2*tcdf(-abs(tmod), dfTotal);
    effect = est;
else
    testname = testname + "F" + blockNote;
    R = covC ./ (unscaled'*unscaled);        % correlation of contrasts
    r = nG - 1;
    F = sum((tmod/R) .* tmod, 2) / r;
    % limma refers the moderated F to the uncapped df (black-box check
    % against speckle), and to chi-square/r when the prior df is infinite.
    if isinf(dfPost)
        pvals = chi2cdf(r*F, r, "upper");
    else
        pvals = fcdf(F, r, dfPost, "upper");
    end
    effect = nan(numel(ctypes), 1);
end
end


function [s2post, dfTotal, dfPost] = i_squeezevar(s2, df)
% Empirical-Bayes variance moderation (Smyth 2004, Section 6.2): fit a
% scaled inverse chi-square prior s0^2, d0 to the sample variances by
% moments of log(s2), then shrink each variance toward s0^2.
n = numel(s2);
m = median(s2);
if m <= 0
    m = 1;
end
s2f = max(s2, 1e-5*m);  % a cell type that is constant in every unit

e = log(s2f) - psi(df/2) + log(df/2);
emean = mean(e);
evar = sum((e - emean).^2)/(n - 1) - psi(1, df/2);
if evar > 0
    d0 = 2*i_trigammainv(evar);
    s02 = exp(emean + psi(d0/2) - log(d0/2));
    s2post = (d0*s02 + df*s2)/(d0 + df);
else
    % No evidence of variance differences between cell types: every cell
    % type gets the pooled variance. (Not exp(emean): that is the moment
    % estimate on the log scale, and differs from the pooled mean.)
    d0 = Inf;
    s2post = repmat(mean(s2), n, 1);
end
% DFPOST is the posterior df. For the t-test it is capped at what pooling
% every cell type would give, since the prior cannot add more than that.
dfPost = df + d0;
dfTotal = min(dfPost, n*df);
end


function y = i_trigammainv(x)
% Solve psi(1, y) == x for y > 0. TRIGAMMA is decreasing from Inf to 0, so
% the root is unique; solve on log(y) to keep the search positive.
g = @(u) psi(1, exp(u)) - x;
lo = log(min(1/x, 1/sqrt(x))) - 2;
hi = log(max(1/x, 1/sqrt(x))) + 2;
y = exp(fzero(g, [lo, hi], optimset(TolX=1e-14)));
end
