function [pvals, info] = i_xctnull(P_s, P_t, li_idx, ri_idx, W12, opts)
%I_XCTNULL  Left-tail null test for candidate L-R pairs, as the reference does.
%   PVALS = TEN.I_XCTNULL(P_S, P_T, LI, RI, W12) scores each candidate pair by
%   the fraction of non-candidate gene pairs whose embedding distance is smaller.
%   A small p-value means the pair sits unusually close in the aligned manifold.
%
%   This reproduces scTenifoldXct's stat.py::null_test. The null is the
%   EXHAUSTIVE set of non-candidate pairs, so the result is deterministic.
%
%   WHY THIS REPLACED A SAMPLED NULL. The previous MATLAB implementation drew
%   max(50000, 100*n_cand) random pairs with randi, which differed from the
%   reference in three ways: it was stochastic, it sampled with replacement, and
%   it did not exclude the candidate pairs from their own null. The first of
%   those was not harmless - Notes.md 6b measured 4 of 31 pairs crossing
%   p = 0.05 on the RNG seed alone, and advised fixing the seed or raising
%   n_null. That advice was treating a symptom: the reference has no such
%   instability because its null is exhaustive.
%
%   NAME-VALUE ARGUMENTS:
%     filterzeros - exclude pairs whose correspondence weight is zero from the
%                   null, as stat.py's filter_zeros=True default does.
%                   "auto" (default) enables it only when W12 is dense.
%
%                   The reference always runs with a dense correspondence (see
%                   Notes.md 6d), where the filter removes almost nothing. Under
%                   this toolbox's sparse "lr" correspondence it would remove
%                   everything: the only nonzero entries ARE the candidates, so
%                   filtering to nonzeros and then excluding candidates leaves an
%                   empty null. "auto" therefore follows the reference where the
%                   comparison is meaningful and disables the filter where it is
%                   degenerate. Pass true or false to override.
%
%   OUTPUTS:
%     pvals - one p-value per candidate pair, in the order of LI_IDX/RI_IDX
%     info  - struct with .n_null, .filterzeros, .n_ties
%
%   TIE HANDLING. scipy's percentileofscore defaults to kind="rank", which
%   averages the ranks of equal values:
%       pct = (n_below + n_at_or_below + (n_at_or_below > n_below)) / (2*n_null)
%   That is reproduced exactly. Embedding distances are continuous, so exact
%   ties are vanishingly rare and pct reduces to n_below/n_null in practice;
%   info.n_ties reports whether any occurred.
%
%   MEMORY. The null requires the full NG-by-NG distance matrix, as the
%   reference's cdist does. That is one quarter the size of the 2*NG-by-2*NG
%   alignment graph the caller already holds when a GRN offset is in use, so it
%   is not the binding constraint on problem size.
%
% see also: TEN.I_XCTCORE, TEN.XCT.XCTMAIN, TEN.I_XCTBLOCK

arguments
    P_s {mustBeNumeric}
    P_t {mustBeNumeric}
    li_idx (:, 1) double
    ri_idx (:, 1) double
    W12 {mustBeNumeric}
    opts.filterzeros = "auto"
end

ng = size(P_s, 1);
if size(P_t, 1) ~= ng
    error("TEN:I_XCTNULL:BadSizes", ...
        "P_S and P_T must have the same number of rows.");
end

% Resolve the filter policy. See the help above for why "auto" exists.
if isstring(opts.filterzeros) || ischar(opts.filterzeros)
    if string(opts.filterzeros) ~= "auto"
        error("TEN:I_XCTNULL:BadFilterZeros", ...
            "filterzeros must be true, false, or ""auto"".");
    end
    filterZeros = nnz(W12) > 0.5*numel(W12);     % dense correspondence
else
    filterZeros = logical(opts.filterzeros);
end

% All-pairs source-to-target distances, as stat.py works from cdist output.
D = pdist2(double(P_s), double(P_t));

candLin = sub2ind([ng, ng], li_idx, ri_idx);
cand_d = D(candLin);

isCand = false(ng, ng);
isCand(candLin) = true;

nullMask = ~isCand;
if filterZeros
    nullMask = nullMask & (W12 ~= 0);
end

null_d = D(nullMask);
nNull = numel(null_d);
if nNull == 0
    error("TEN:I_XCTNULL:EmptyNull", ...
        "The null distribution is empty. This happens when filterzeros is " + ...
        "true and the correspondence block is sparse, because its only " + ...
        "nonzero entries are the candidate pairs themselves. Pass " + ...
        "filterzeros=false, or use a dense correspondence (w12mode=""outer"").");
end

pvals = i_percentileofscore(null_d, cand_d, nNull);

info = struct("n_null", nNull, "filterzeros", filterZeros, ...
    "n_ties", sum(ismember(cand_d, null_d)));

end


%% ---- scipy.stats.percentileofscore, kind="rank" ----
function pct = i_percentileofscore(null_d, cand_d, nNull)

nullSorted = sort(null_d(:));
[candSorted, order] = sort(cand_d(:));

% Counts of null values strictly below each candidate. histcounts bins are
% half-open [e(k), e(k+1)), so the running total up to bin k is exactly
% #{null < candSorted(k)}.
edges = [-Inf; candSorted; Inf];
below = cumsum(histcounts(nullSorted, edges));
below = below(1:numel(candSorted))';

% Exact ties. Vanishingly rare for continuous distances, so the membership
% test short-circuits in the normal case.
ties = zeros(numel(candSorted), 1);
isTied = ismember(candSorted, nullSorted);
for k = find(isTied)'
    ties(k) = sum(nullSorted == candSorted(k));
end
atOrBelow = below + ties;

pctSorted = (below + atOrBelow + (atOrBelow > below))./(2*nNull);

pct = zeros(numel(candSorted), 1);
pct(order) = pctSorted;

end
