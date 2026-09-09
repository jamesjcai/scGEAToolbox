function [pval, sigma] = e_deviationpvalue(d)
%E_DEVIATIONPVALUE  P-value for a non-negative deviation statistic.
%
%   pval = PKG.E_DEVIATIONPVALUE(d) returns, for each entry of the
%   non-negative statistic D, the probability of seeing a deviation at
%   least that large under an empirical null fitted to D itself.
%
%   [pval, sigma] = ... also returns the fitted null scale.
%
%   THE TWO THINGS THIS GETS RIGHT. Both were wrong in the callers, and
%   both push the same way, so together they inflated the number of hits
%   several-fold.
%
%   1. D is a magnitude -- a distance, a norm, an absolute difference. Its
%      null is the folded half of a symmetric distribution centred at zero,
%      so the p-value has to account for both tails of that distribution:
%      P(|Z| > d) = 2*P(Z > d). SC_SPLINEFIT fitted a symmetric null by
%      mirroring its sample and then read one tail of it, which halved
%      every p-value it reported.
%
%   2. The null scale comes from the median, not from the standard
%      deviation of a truncated bulk. Truncating a sample and then taking
%      its SD underestimates the scale badly: on 20000 draws from an exact
%      half-normal null with sigma 1.400, the sub-90th-percentile SD
%      returned 1.108. Too small a scale is anti-conservative on top of
%      the missing factor of two. For d ~ |N(0, sigma)| the median is
%      norminv(0.75)*sigma, and a median is already immune to the upper
%      tail that carries the signal, so nothing has to be truncated away.
%
%   Measured on that half-normal null, P(p < 0.05) was 0.193 as the
%   callers had it, 0.122 with the factor of two alone, and 0.052 with
%   both -- against a nominal 0.05.
%
%   THE FLOOR AT 1/N. A null fitted from N observations cannot resolve a
%   p-value far below 1/N, so the result is floored there. Without it the
%   normal tail gets extrapolated to absurdity: the far tail of these
%   statistics is heavier than normal -- on a null split the 99.9th
%   percentile sat at 4.6 sigma where a half-normal predicts 3.3 -- and a
%   handful of genes came back at 1e-42 from a 2000-point fit. Those were
%   not close calls; they were enough to put 11 genes past BH q<=0.05 on
%   one homogeneous population split at random, where a calibrated test
%   should call nothing.
%
%   The floor costs no power, because the genes that reach it reach it
%   together and BH clears a block of ties easily. Measured on a contrast
%   with 400 planted genes among 2000: floored, 406 called at precision
%   0.96 and recall 0.97; unfloored, the same 406 at the same precision
%   and recall -- and 11 false calls on the null split instead of none.
%   Ranking among the tied top hits is not lost either: it lives in the
%   statistic itself, which is what SC_SPLINEFIT and SC_DVG sort on.
%
%   ON WHAT THIS CAN AND CANNOT SAY. The null is fitted to the observed
%   statistic, so it assumes most entries are null -- the usual empirical
%   null. That holds for deviation-from-trend and differential-variability
%   statistics, where a handful of genes depart and the rest do not. It
%   does not hold if the effect is global: if every gene were more variable
%   in one group, the median would move with them and nothing would be
%   called. A statistic of that shape needs a permutation null instead.
%   Fitting the tail rather than the body is not an option here for the
%   same reason: on a real contrast the tail IS the signal.
%
%   INPUT:
%     d - vector of non-negative deviations. Zeros are kept and get
%         pval = 1; NaNs propagate.
%
%   OUTPUTS:
%     pval  - same size as D, in [1/numel(D), 1].
%     sigma - the fitted scale, or NaN when D carries no usable spread.
%
%   See also SC_SPLINEFIT, SC_DVG, SC_HVG.

arguments
    d {mustBeNumeric}
end

pval = ones(size(d));
finite = isfinite(d);
sigma = NaN;
if ~any(finite)
    pval(~finite) = NaN;
    return;
end

% norminv(0.75) is the median of |N(0,1)|, so median(d)/norminv(0.75)
% estimates the scale of the symmetric distribution d is the folded half of.
sigma = median(double(d(finite)))/norminv(0.75);

% A median of zero means over half the entries are exactly zero, which
% happens when a caller has flattened the uninteresting ones. Fall back to
% the median over what is left rather than reporting a degenerate scale.
if sigma <= 0
    positive = finite & d > 0;
    if any(positive)
        sigma = median(double(d(positive)))/norminv(0.75);
    end
end

if ~isfinite(sigma) || sigma <= 0
    pval(~finite) = NaN;
    return;
end

pval(finite) = min(1, 2*normcdf(double(d(finite)), 0, sigma, "upper"));
pval(finite) = max(pval(finite), 1/nnz(finite));
pval(~finite) = NaN;

end
