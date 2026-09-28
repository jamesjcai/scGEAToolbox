function [T] = e_findallmarkers(X, g, c, cL, logfc, minpct, showwaitbar, ...
                maxnummarkers, opts)
% E_FINDALLMARKERS  Markers of every cluster against the rest, as Seurat's FindAllMarkers.
%   T = pkg.e_findallmarkers(X, g, c) tests each group in C against all other
%   cells with a Wilcoxon rank-sum test on log1p-normalised counts.
%
%   Positional options, [] for the default:
%     cL             group names, when C is already integer-coded
%     logfc          |avg_log2FC| >= logfc              (logfc.threshold, 0.1)
%     minpct         max(pct_1, pct_2) >= minpct        (min.pct, 0.01)
%     showwaitbar    false
%     maxnummarkers  rows kept per group after ranking  (Inf)
%
%   Name-value options:
%     OnlyPos       keep only avg_log2FC > 0 rows       (only.pos, false)
%     PAdjust       "bonferroni" over all genes, as Seurat, or "bh"
%     ReturnThresh  significance cut-off                (return.thresh, 0.01)
%     ThresholdOn   column ReturnThresh applies to: "p_val", as Seurat, or
%                   "p_val_adj"
%
%   Every default is Seurat v5's. The pre-v5 behaviour of this function was
%   logfc 0.5, minpct 0.1, OnlyPos=true, PAdjust="bh", ReturnThresh=0.05,
%   ThresholdOn="p_val_adj".
%
%   https://satijalab.org/seurat/reference/findallmarkers
arguments
    X
    g
    c
    cL = []
    logfc = []
    minpct = []
    showwaitbar = false
    maxnummarkers = []
    opts.OnlyPos (1,1) logical = false
    opts.PAdjust (1,1) string {mustBeMember(opts.PAdjust, ["bonferroni", "bh"])} = "bonferroni"
    opts.ReturnThresh (1,1) double {mustBeBetween(opts.ReturnThresh, 0, 1)} = 0.01
    opts.ThresholdOn (1,1) string {mustBeMember(opts.ThresholdOn, ["p_val", "p_val_adj"])} = "p_val"
end
% Callers pass [] to mean "default", which an arguments-block default does
% not cover.
if isempty(logfc), logfc = 0.1; end
if isempty(minpct), minpct = 0.01; end
if isempty(maxnummarkers), maxnummarkers = Inf; end

if isempty(cL)
    [c, cL] = findgroups(string(c));
end
% Kept sparse. This used to FULL() the whole matrix and then copy
% X(:, c == kc) and X(:, c ~= kc) for every cluster before a per-gene
% RANKSUM loop: 37.9 s and a dense copy of the matrix at 5000 genes x
% 10000 cells x 10 clusters, over 20 GB at 20000 x 50000.
X = log1p(sc_norm(X));
if showwaitbar
    fw = gui.myWaitbar([]);
end
mC = max(c);
c = c(:);
N = size(X, 2);

% Every cluster against the rest from one ranking of each gene; the same
% p-values RANKSUM gives. See PKG.E_RANKSUMROWS.
P = pkg.e_ranksumrows(X, c);

% Means on the normalised count scale and detection rates, per cluster and
% for the rest, through a sparse cells-by-clusters membership matrix.
% EXPM1 inverts the LOG1P above exactly and keeps X sparse.
member = sparse(1:N, c, 1, N, mC);
nIn = full(sum(member, 1));
nOut = N - nIn;
E = expm1(X);
sumIn = full(E*member);
avgIn = sumIn./nIn;
avgOut = (full(sum(E, 2)) - sumIn)./nOut;
detIn = full(double(X > 0)*member);
pctIn = detIn./nIn;
pctOut = (full(sum(X > 0, 2)) - detIn)./nOut;

Tcell = cell(mC, 1);
for kc = 1:mC
    if showwaitbar
        if kc ~= mC
            gui.myWaitbar([], fw, [], [], sprintf('Processing %s', cL{kc}), kc/mC);
        else
            gui.myWaitbar([], fw, [], [], sprintf('Processing %s', cL{kc}), (kc-1)/mC);
        end
    end
    [t] = in_findmarkers(P(:, kc), avgIn(:, kc), avgOut(:, kc), ...
        pctIn(:, kc), pctOut(:, kc), cL(kc));
    % IN_FINDMARKERS returns its rows already ranked. It has to: truncating
    % first and sorting afterwards -- which is what this did -- keeps the
    % first MAXNUMMARKERS genes in GENELIST order and throws the real
    % markers away. On a 3000-gene fixture with the 100 strongest markers
    % placed late in the list and 300 weak ones early, the returned "top
    % 100" contained none of the strong ones and all 100 of the weak.
    Tcell{kc} = t(1:min([maxnummarkers, size(t,1)]), :);
end
T = vertcat(Tcell{:});
[~, idx] = sort(T.p_val_adj);
T = T(idx, :);
[~, idx] = natsort(T.grp);
T = T(idx, :);
if showwaitbar
    gui.myWaitbar([], fw);
end


function [t] = in_findmarkers(p_val, avg_1, avg_2, pct_1, pct_2, ctxt)
ng = numel(p_val);
% AVG_1 and AVG_2 are means on the normalised count scale, not of the
% log1p values. Taking log2 of a ratio of log-means -- which this once
% computed -- is not a fold change on any scale: a marker at 200 against 70
% counts has a true log2FC of 1.50 and came out as 0.315, below the default
% LOGFC threshold, so the strongest markers were the ones being dropped.
% Pseudocount and form follow SC_DEG so that the two functions report the
% same number for the same contrast.
avg_log2FC = log2(avg_1 + 1) - log2(avg_2 + 1);
if opts.PAdjust == "bonferroni"
    % Seurat's p.adjust(..., n = nrow(object)): the family is every gene,
    % tested or not. A NaN p-value stays NaN and fails every threshold.
    p_val_adj = min(p_val*ng, 1);
else
    % PKG.E_FDR replaces the two-branch block that used to sit here. Its
    % MAFDR branch and its PKG.E_FDR_BH branch are the same algorithm on
    % clean input -- they agree to 2e-16 -- but they disagree whenever the
    % p-values carry NaN, which is what RANKSUM returns for a gene with no
    % counts in either group and what a per-cell-type run produces in bulk.
    % One drops those from the family, the other counts them, so the same
    % command gave a different answer depending on whether the
    % Bioinformatics Toolbox was installed.
    p_val_adj = pkg.e_fdr(p_val);
end
grp = repmat(ctxt, ng, 1);
t = table(grp, g, p_val, avg_log2FC, avg_1, avg_2, ...
    pct_1, pct_2, p_val_adj);
if opts.OnlyPos
    isFoldChangeKept = (avg_log2FC >= logfc) & (avg_log2FC > 0);
else
    isFoldChangeKept = abs(avg_log2FC) >= logfc;
end
t = t(t.(opts.ThresholdOn) < opts.ReturnThresh & isFoldChangeKept & ...
    max(pct_1, pct_2) >= minpct, :);
% Rank here, before the caller truncates to MAXNUMMARKERS. Ties on the
% adjusted p-value are common -- the rank-sum test has a floor set by the
% group sizes -- so the fold change breaks them, as in SC_DEGMAST.
t = sortrows(t, ["p_val_adj", "avg_log2FC"], ["ascend", "descend"]);
end

end
