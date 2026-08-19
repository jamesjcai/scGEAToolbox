function [pvals, qvals, FC] = i_xctchi2(stat, opts)
%I_XCTCHI2  Chi-square enrichment test on aligned distances, as the reference.
%   [PVALS, QVALS, FC] = TEN.I_XCTCHI2(STAT) scores each pair by referring its
%   statistic to a chi-square distribution after scaling by the mean statistic:
%
%       FC = dof * STAT / mean(STAT)
%       p  = chi2cdf(FC, dof)        for tail="left"
%       p  = 1 - chi2cdf(FC, dof)    for tail="right"
%
%   This reproduces stat.py's chi2_test (tail="left", STAT = dist^2) and
%   chi2_diff_test (tail="right", STAT = diff2 = (dist_A - dist_B)^2).
%
%   STAT MUST COVER ALL GENE PAIRS, not just the candidates. The scaling is by
%   `mean(STAT)`, and the reference computes that over its full all-pairs table;
%   restricting to candidates would change the denominator and hence every
%   p-value. The same applies to the FDR step, whose correction factor is the
%   number of pairs tested.
%
%   NAME-VALUE ARGUMENTS:
%     dof  - chi-square degrees of freedom (default 1, as the reference)
%     tail - "right" (default) for a difference statistic, where large values
%            are significant; "left" for a distance statistic, where small
%            values are significant
%     fdr  - apply Benjamini-Hochberg correction (default true, as the
%            reference's FDR=True). QVALS is empty when false.
%
%   OUTPUTS:
%     pvals - raw p-value per pair
%     qvals - BH-adjusted q-value per pair, or [] when fdr=false
%     FC    - the scaled statistic actually referred to the chi-square
%
%   THE FDR IS OVER EVERY PAIR, NOT THE CANDIDATES. stat.py calls
%   multipletests on the full p-vector and only afterwards restricts to
%   candidates, so the correction factor is the total number of gene pairs. That
%   is a strong correction - millions of hypotheses - and it is deliberate on
%   the reference's part. Filtering to candidates first would be far more
%   permissive and would not reproduce its numbers.
%
%   BH is implemented directly rather than via mafdr, to avoid a Bioinformatics
%   Toolbox dependency. It matches statsmodels' fdr_bh, including the
%   monotonicity enforcement that makes q non-decreasing in p.
%
% see also: TEN.I_XCTNULL, TEN.SCTENIFOLDXCT2, TEN.XCT.XCTMAIN2

arguments
    stat {mustBeNumeric}
    opts.dof (1, 1) double {mustBePositive} = 1
    opts.tail (1, 1) string {mustBeMember(opts.tail, ["left", "right"])} = "right"
    opts.fdr (1, 1) logical = true
end

stat = double(stat(:));
m = numel(stat);
if m == 0
    error("TEN:I_XCTCHI2:EmptyStat", "STAT is empty.");
end

statMean = mean(stat);
if statMean <= 0
    error("TEN:I_XCTCHI2:DegenerateStat", ...
        "mean(STAT) is %g; the chi-square scaling is undefined. All pairs " + ...
        "having an identical statistic usually means the embedding collapsed.", ...
        statMean);
end

FC = opts.dof.*stat./statMean;

if opts.tail == "left"
    pvals = chi2cdf(FC, opts.dof);
else
    pvals = 1 - chi2cdf(FC, opts.dof);
end

if opts.fdr
    qvals = i_bhfdr(pvals);
else
    qvals = [];
end

end


%% ---- Benjamini-Hochberg, matching statsmodels' fdr_bh ----
function q = i_bhfdr(p)
m = numel(p);
[pSorted, order] = sort(p(:));

% q_(i) = min over k >= i of (m/k)*p_(k), enforced by a reverse running
% minimum, then clipped at 1. cummin(flip(...)) is the vectorised form of
% statsmodels' np.minimum.accumulate on the reversed array.
raw = pSorted.*(m./(1:m)');
qSorted = flip(cummin(flip(raw)));
qSorted = min(qSorted, 1);

q = zeros(m, 1);
q(order) = qSorted;
q = reshape(q, size(p));
end
