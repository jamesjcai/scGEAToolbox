function [T, info] = sc_gsettest(stats, genelist, setmatrx, setnames, setgenes, opts)
%SC_GSETTEST  Competitive gene-set tests for a ranked gene list.
%
%   T = SC_GSETTEST(stats, genelist, setmatrx, setnames, setgenes) tests
%   every gene set against the ranked list defined by STATS using a
%   competitive null hypothesis: "the genes in set S are no more extreme
%   than the genes outside S". This is the null behind GSEA, CAMERA and
%   over-representation analysis, and it differs from the self-contained
%   null ("no gene in S is differentially expressed") tested by ROAST.
%
%   All methods share one backend: a sparse nSets-by-nGenes membership
%   matrix built once over the measured gene universe. Every closed-form
%   method then reduces to a sparse matrix-vector product, so a full MSigDB
%   collection is screened in well under a second with no permutation loop
%   and no Python/R dependency.
%
%   USAGE:
%     [setmatrx, setnames, setgenes] = pkg.e_getgenesets('MSIGDB', 'human');
%     Tdeg = sc_deg(X1, X2, genelist);
%     T = sc_gsettest(Tdeg.tScore, Tdeg.gene, setmatrx, setnames, setgenes);
%     T = sc_gsettest(s, g, setmatrx, setnames, setgenes, Method="camera");
%     T = sc_gsettest(s, g, {["A" "B" "C"], ["D" "E"]}, ["s1" "s2"], []);
%
%   INPUTS:
%     stats    - G-by-1 gene-level statistic defining the ranking, e.g. a
%                moderated t-statistic, log2 fold change, or signed
%                -log10(p). Larger = "more up". NaNs are dropped.
%     genelist - G-by-1 gene symbols, aligned with STATS.
%     setmatrx - nSets-by-nSetGenes logical or numeric membership matrix, as
%                returned by PKG.E_GETGENESETS. Numeric values are treated
%                as gene weights by Method="ulm" and as plain membership by
%                all other methods. May instead be an nSets-by-1 cell array
%                of gene-name lists, in which case SETGENES is ignored.
%     setnames - nSets-by-1 set names ([] to auto-generate).
%     setgenes - nSetGenes-by-1 gene symbols for the columns of SETMATRX.
%
%   NAME-VALUE ARGUMENTS:
%     Method       - "ranksum" (default) Wilcoxon rank-sum of set vs rest,
%                      tie-corrected normal approximation. limma's
%                      geneSetTest/wilcoxGST. Usually more powerful than the
%                      weighted KS statistic of GSEA and needs no permutation.
%                    "camera" Two-sample t on STATS with the variance
%                      inflated by VIF = 1 + (m-1)*rho to absorb inter-gene
%                      correlation (Wu & Smyth 2012). The honest competitive
%                      test - prefer it whenever co-expression matters.
%                    "cameraRank" Rank-sum variant of CAMERA (the
%                      use.ranks=TRUE form); robust to non-normal STATS.
%                    "cerno" Fisher's method on ranks from the top,
%                      -2*sum(log(r/n)) ~ chi2 with 2m df (tmod). Fast, and
%                      sensitive to enrichment concentrated at the very top.
%                    "page" Parametric Analysis of Gene set Enrichment
%                      z-score (Kim & Volsky 2005).
%                    "ulm" Univariate linear model: regress STATS on set
%                      membership and report the t-value of the slope. For
%                      binary sets this equals "camera"; its value is that
%                      it accepts signed/weighted membership, e.g. DoRothEA
%                      regulons carrying a mode of regulation.
%                    "gsea" Weighted Kolmogorov-Smirnov enrichment score
%                      with gene-label permutation. The slow path, provided
%                      for comparability with published GSEA results.
%     Direction    - "both" (default) two-sided, "up", "down", or "mixed"
%                    (test |STATS|, i.e. "the set is perturbed in either
%                    direction"). Reporting "both" alongside "mixed" is the
%                    convention romer follows and is usually informative.
%     MinSize      - minimum genes of a set present in GENELIST (default 5).
%     MaxSize      - maximum (default Inf; 500 is a common choice).
%     InterGeneCor - rho for the CAMERA variance inflation factor, a scalar
%                    or one value per set (default 0.01, limma's default).
%                    Set 0 to recover the uncorrected test. Used only by
%                    "camera", "cameraRank", "page" and "ulm".
%     Weight       - exponent on |STATS| for Method="gsea" (default 1;
%                    0 gives the unweighted, classic KS statistic).
%     NumPerm      - permutations for Method="gsea" (default 1000).
%     Sort         - sort rows by ascending PValue (default true).
%     Verbose      - print a short summary (default true).
%
%   OUTPUTS:
%     T    - table with one row per tested set:
%              SetName  set name
%              SetSize  genes of the set found in GENELIST
%              Stat     the method's native statistic, signed so that
%                       positive means the set sits toward the top
%              AUC      Mann-Whitney effect size, comparable across every
%                       method. 0.5 = no shift, >0.5 = set ranks higher
%              PValue   raw p-value
%              FDR      Benjamini-Hochberg adjusted p-value
%            Method="gsea" adds ES and NES columns.
%     info - struct holding the membership matrix, gene universe, working
%            statistic, and the options actually used.
%
%   NOTE ON CORRELATION. Every competitive test permutes genes, which
%   assumes gene independence; co-expression breaks that and inflates the
%   FDR. "camera" is the only method here that corrects for it. For
%   single-cell data the problem is compounded by pseudo-replication:
%   thousands of cells produce enormous gene-level statistics and every set
%   looks significant. Compute STATS on sample-level pseudobulk rather than
%   per cell before running any of these tests.
%
% REF: Goeman & Buhlmann (2007) Bioinformatics 23:980 (competitive vs
%      self-contained); Wu & Smyth (2012) NAR 40:e133 (CAMERA);
%      Kim & Volsky (2005) BMC Bioinformatics 6:144 (PAGE);
%      Zyla et al. (2019) Bioinformatics 35:5146 (CERNO);
%      Subramanian et al. (2005) PNAS 102:15545 (GSEA).
%
% See also: SC_DEG, SC_DVG, SC_HVG, PKG.E_GETGENESETS, PKG.E_ULM

arguments
    stats {mustBeNumeric, mustBeVector}
    genelist (:, 1) string
    setmatrx
    setnames = []
    setgenes = []
    opts.Method (1, 1) string {mustBeMember(opts.Method, ["ranksum", ...
        "camera", "cameraRank", "cerno", "page", "ulm", "gsea"])} = "ranksum"
    opts.Direction (1, 1) string {mustBeMember(opts.Direction, ...
        ["both", "up", "down", "mixed"])} = "both"
    opts.MinSize (1, 1) double = 5
    opts.MaxSize (1, 1) double = Inf
    opts.InterGeneCor (:, 1) double = 0.01
    opts.Weight (1, 1) double = 1
    opts.NumPerm (1, 1) double = 1000
    opts.Sort (1, 1) logical = true
    opts.Verbose (1, 1) logical = true
end

stats = double(stats(:));
if numel(stats) ~= numel(genelist)
    error("sc_gsettest:SizeMismatch", ...
        "STATS (%d) and GENELIST (%d) must have the same length.", ...
        numel(stats), numel(genelist));
end

[setmatrx, setnames, setgenes] = i_normalizesets(setmatrx, setnames, setgenes);

% Drop unusable genes, then collapse duplicate symbols to their most extreme
% statistic. The universe is every measured gene, including genes in no set:
% those still form part of the "rest" the set is competed against.
keep = isfinite(stats) & strlength(genelist) > 0;
stats = stats(keep);
genelist = genelist(keep);

[ug, ~, ic] = unique(upper(genelist));
if numel(ug) < numel(genelist)
    ndup = numel(genelist) - numel(ug);
    grp = accumarray(ic, (1:numel(stats))', [], @(v) {v});
    sel = cellfun(@(v) i_argmaxabs(stats, v), grp);
    warning("sc_gsettest:DuplicateGenes", ...
        "%d duplicate gene symbols collapsed to their most extreme statistic.", ...
        ndup);
    stats = stats(sel);
    genelist = genelist(sel);
end

n = numel(stats);
if n < 10
    error("sc_gsettest:TooFewGenes", ...
        "Only %d usable genes; a competitive test is meaningless.", n);
end

% Shared backend: one sparse membership matrix over the gene universe.
if ~issparse(setmatrx)
    setmatrx = sparse(double(setmatrx));
end
[tf, loc] = ismember(upper(genelist), upper(setgenes));
M = sparse(size(setmatrx, 1), n);
M(:, tf) = setmatrx(:, loc(tf));

msize = full(sum(M ~= 0, 2));
ok = msize >= opts.MinSize & msize <= min(opts.MaxSize, n - 1);
if ~any(ok)
    error("sc_gsettest:NoSets", ...
        "No gene set has between %g and %g of its genes in GENELIST " + ...
        "(largest overlap was %d). Check that GENELIST and SETGENES use " + ...
        "the same symbol namespace.", opts.MinSize, opts.MaxSize, max(msize));
end
W = M(ok, :);                   % weighted form, kept for "ulm"
M = spones(W);                  % plain membership for everything else
names = setnames(ok);
m = full(sum(M, 2));
nSets = numel(m);

switch opts.Direction
    case "both"
        s = stats;
        tail = "two";
    case "up"
        s = stats;
        tail = "one";
    case "down"
        s = -stats;
        tail = "one";
    case "mixed"
        s = abs(stats);
        tail = "one";
end

rho = opts.InterGeneCor;
if isscalar(rho)
    rho = repmat(rho, nSets, 1);
elseif numel(rho) == numel(ok)
    rho = rho(ok);
elseif numel(rho) ~= nSets
    error("sc_gsettest:BadRho", ...
        "InterGeneCor must be a scalar or have one value per gene set.");
end
vif = max(1 + (m - 1).*rho, 1e-12);

% AUC is reported for every method so effect sizes stay comparable.
r = tiedrank(s);
R = full(M*r);
auc = (R - m.*(m + 1)/2) ./ (m.*(n - m));

extra = table();

switch opts.Method
    case {"ranksum", "cameraRank"}
        [~, ~, ic2] = unique(s);
        tcnt = accumarray(ic2, 1);
        tieadj = sum(tcnt.^3 - tcnt);
        muR = m*(n + 1)/2;
        varR = m.*(n - m)/12 .* ((n + 1) - tieadj/(n*(n - 1)));
        if opts.Method == "cameraRank"
            varR = varR.*vif;
        end
        d = R - muR;
        stat = sign(d).*max(abs(d) - 0.5, 0) ./ sqrt(varR);
        p = i_normtail(stat, tail);

    case "camera"
        [stat, df] = i_twosamplet(M, s, m, n);
        stat = stat ./ sqrt(vif);
        p = i_ttail(stat, df, tail);

    case "page"
        xin = full(M*s) ./ m;
        stat = (xin - mean(s)).*sqrt(m) ./ (std(s).*sqrt(vif));
        p = i_normtail(stat, tail);

    case "ulm"
        % Regress s on (possibly signed or weighted) membership and report
        % the t-value of the slope, all sets at once.
        wsum = full(sum(W, 2));
        sxx = full(sum(W.^2, 2)) - wsum.^2/n;
        sxx(sxx <= 0) = NaN;
        syy = sum((s - mean(s)).^2);
        sxy = full(W*s) - wsum*mean(s);
        b1 = sxy ./ sxx;
        sse = max(syy - b1.*sxy, 0);
        df = n - 2;
        stat = b1 ./ sqrt((sse/df) ./ sxx .* vif);
        p = i_ttail(stat, df, tail);

    case "cerno"
        [stat, p] = i_cerno(M, s, m, n, tail);

    case "gsea"
        [stat, p, es] = i_gsea(M, s, m, n, opts);
        extra = table(es, stat, VariableNames=["ES", "NES"]);
end

[~, ~, ~, fdr] = pkg.e_fdr_bh(p, 0.05, 'pdep', 'no');

T = table(string(names(:)), m, stat, auc, p, fdr(:), VariableNames= ...
    ["SetName", "SetSize", "Stat", "AUC", "PValue", "FDR"]);
if ~isempty(extra)
    T = [T, extra];
end
if opts.Sort
    T = sortrows(T, "PValue");
end

if opts.Verbose
    fprintf("sc_gsettest: method=%s, direction=%s\n", opts.Method, opts.Direction);
    fprintf("  universe %d genes; %d of %d sets tested (size %d-%d)\n", ...
        n, nSets, numel(ok), min(m), max(m));
    if ismember(opts.Method, ["camera", "cameraRank", "page", "ulm"])
        fprintf("  inter-gene correlation rho=%g (VIF %.2f-%.2f)\n", ...
            rho(1), min(vif), max(vif));
    end
    fprintf("  %d sets at FDR<0.05, %d at FDR<0.25\n", ...
        sum(T.FDR < 0.05), sum(T.FDR < 0.25));
end

if nargout > 1
    info = struct(M=M, W=W, Universe=genelist, Stat=s, SetIndex=find(ok), ...
        Rho=rho, VIF=vif, Options=opts);
end
end

% =========================================================================
function idx = i_argmaxabs(stats, v)
% Index of the most extreme statistic within one group of duplicate symbols.
[~, k] = max(abs(stats(v)));
idx = v(k);
end

% =========================================================================
function [setmatrx, setnames, setgenes] = i_normalizesets(setmatrx, setnames, setgenes)
% Accept either the PKG.E_GETGENESETS triple or a cell array of gene lists.
if iscell(setmatrx) && ~isempty(setmatrx) && ~isnumeric(setmatrx{1})
    lists = cellfun(@(v) string(v(:)), setmatrx(:), UniformOutput=false);
    setgenes = unique(vertcat(lists{:}));
    setgenes = setgenes(strlength(setgenes) > 0);
    mat = false(numel(lists), numel(setgenes));
    for k = 1:numel(lists)
        mat(k, :) = ismember(setgenes, lists{k});
    end
    setmatrx = mat;
end
if isempty(setgenes)
    error("sc_gsettest:NoSetGenes", ...
        "SETGENES is empty. Pass the gene symbols for the columns of " + ...
        "SETMATRX, or pass SETMATRX as a cell array of gene-name lists.");
end
setgenes = string(setgenes(:));
if isempty(setnames)
    setnames = "Set" + string((1:size(setmatrx, 1))');
end
setnames = string(setnames(:));
if numel(setnames) ~= size(setmatrx, 1)
    error("sc_gsettest:SizeMismatch", ...
        "SETNAMES (%d) must have one entry per row of SETMATRX (%d).", ...
        numel(setnames), size(setmatrx, 1));
end
if numel(setgenes) ~= size(setmatrx, 2)
    error("sc_gsettest:SizeMismatch", ...
        "SETGENES (%d) must have one entry per column of SETMATRX (%d).", ...
        numel(setgenes), size(setmatrx, 2));
end
end

% =========================================================================
function [t, df] = i_twosamplet(M, s, m, n)
% Pooled-variance two-sample t of set genes against the remaining genes,
% for all sets at once via two sparse matrix-vector products.
sum1 = sum(s);
sum2 = sum(s.^2);
Ms = full(M*s);
xin = Ms ./ m;
xout = (sum1 - Ms) ./ (n - m);
sswithin = sum2 - m.*xin.^2 - (n - m).*xout.^2;
df = n - 2;
sp2 = max(sswithin, 0)/df;
t = (xin - xout) ./ sqrt(sp2 .* (1./m + 1./(n - m)));
end

% =========================================================================
function [stat, p] = i_cerno(M, s, m, n, tail)
% Fisher's method on ranks measured from the top of the list.
rtop = tiedrank(-s);
lup = -2*full(M*log(rtop/n));
pup = chi2cdf(lup, 2*m, "upper");
if tail == "one"
    stat = lup;
    p = pup;
    return;
end
% CERNO is one-sided by construction; the two-sided form takes the better
% of the two tails and pays the usual factor of two for looking at both.
ldn = -2*full(M*log(tiedrank(s)/n));
pdn = chi2cdf(ldn, 2*m, "upper");
isup = pup <= pdn;
stat = lup.*isup - ldn.*~isup;
p = min(1, 2*min(pup, pdn));
end

% =========================================================================
function [nes, p, es] = i_gsea(M, s, m, n, opts)
% Preranked GSEA with gene-label permutation. Only the running sum at hit
% positions can be extremal, so each enrichment score costs O(m log m)
% rather than O(n) -- without that the permutation loop is unusable.
[ssort, ord] = sort(s, "descend");
w = abs(ssort).^opts.Weight;
Msort = M(:, ord);
nSets = size(M, 1);
nperm = opts.NumPerm;

es = zeros(nSets, 1);
nes = zeros(nSets, 1);
p = ones(nSets, 1);

for k = 1:nSets
    mk = m(k);
    es(k) = i_es(w, find(Msort(k, :))', n, mk);

    esp = zeros(nperm, 1);
    for b = 1:nperm
        esp(b) = i_es(w, sort(randperm(n, mk))', n, mk);
    end

    if es(k) >= 0
        pool = esp(esp > 0);
        nhit = sum(esp >= es(k));
    else
        pool = esp(esp < 0);
        nhit = sum(esp <= es(k));
    end
    if ~isempty(pool)
        nes(k) = es(k)/abs(mean(pool));
    end
    p(k) = (nhit + 1)/(nperm + 1);
end

% "up", "down" and "mixed" all reduce to top-of-list enrichment once S has
% been transformed, so scores pointing the other way are not evidence.
if opts.Direction ~= "both"
    p(es < 0) = 1;
end
end

% =========================================================================
function es = i_es(w, q, n, m)
% Weighted KS enrichment score. Q holds the ascending positions of the M set
% genes within the descending-statistic ordering.
nr = sum(w(q));
if nr <= 0 || m >= n
    es = 0;
    return;
end
cumhit = cumsum(w(q))/nr;
misses = (q - (1:m)')/(n - m);
espos = max([0; cumhit - misses]);                  % maxima land on a hit
esneg = min([0; [0; cumhit(1:m-1)] - misses]);      % minima just before one
if espos > -esneg
    es = espos;
else
    es = esneg;
end
end

% =========================================================================
function p = i_normtail(z, tail)
if tail == "two"
    p = 2*normcdf(-abs(z));
else
    p = normcdf(z, "upper");
end
end

% =========================================================================
function p = i_ttail(t, df, tail)
if tail == "two"
    p = 2*tcdf(-abs(t), df);
else
    p = tcdf(t, df, "upper");
end
end
