function [thres, info] = e_bimodalthres(scores, opts)
%E_BIMODALTHRES  Split a bimodal score distribution at its valley.
%
%   thres = PKG.E_BIMODALTHRES(scores) returns the value that separates the
%   two modes of SCORES, or [] when the distribution is not convincingly
%   bimodal.
%
%   The distribution is smoothed with a Gaussian kernel, its local extrema
%   are located, and the threshold is placed at the lowest point of the
%   density between the two tallest peaks. Unlike a fixed quantile or a
%   two-component Gaussian fit, this makes no assumption about the shape or
%   the relative size of the two populations, and it declines to answer
%   when there is only one population to begin with -- which is the case
%   that matters, because forcing a split on a unimodal score is how
%   spurious "positive" cell populations get invented.
%
%   NAME-VALUE:
%     NumPoints    - grid points for the density estimate (default 512).
%     Bandwidth    - kernel bandwidth. Default is Silverman's rule of thumb,
%                    0.9*min(sd, IQR/1.34)*n^(-1/5).
%     MergeTol     - extrema whose densities differ by less than this
%                    fraction of the tallest peak are treated as one flat
%                    stretch and discarded (default 0.001).
%     ProminenceTol- require the second peak's prominence to be at least
%                    this fraction of the tallest peak's (default 0.05).
%     MinDipDepth  - require the valley to sit at least this far below the
%                    lower of the two peaks, as a fraction of that peak's
%                    height (default 0.05). Set both tolerances to 0 to
%                    accept any pair of extrema.
%
%   OUTPUTS:
%     thres - the cut point, or [] if no split is supported.
%     info  - struct with the density grid (x, y), the extrema found, and
%             the bandwidth used. Handy for plotting the decision.
%
%   DEFAULTS. scCancer uses ProminenceTol = 0.01 and has no dip criterion,
%   which splits unimodal data far too readily. Measured over 300 draws of
%   400 points each, splitting unimodal samples it should have refused:
%
%                                gauss   expon   lognorm | clean  overlap
%     scCancer (0.01, none)        11%     39%       43% |  100%      96%
%     defaults here (0.05, 0.05)  1.3%    0.3%        0% |  100%      82%
%
%   where "clean" is N(0,1) vs N(6,1) and "overlap" is N(0,1) vs N(3,1),
%   both 300:200. The trade is a modest loss of sensitivity on badly
%   overlapping mixtures for a sixty-fold drop in invented splits. Pass
%   ProminenceTol=0.01, MinDipDepth=0 to reproduce scCancer exactly.
%
% REF: Guo et al. (2021) Brief Bioinform 22:bbaa127 (scCancer).
%
% See also: SC_MALIGNSCORE, KSDENSITY.

arguments
    scores {mustBeNumeric, mustBeVector}
    opts.NumPoints (1, 1) double {mustBePositive} = 512
    opts.Bandwidth (1, 1) double {mustBeNonnegative} = 0
    opts.MergeTol (1, 1) double {mustBeNonnegative} = 0.001
    opts.ProminenceTol (1, 1) double {mustBeNonnegative} = 0.05
    opts.MinDipDepth (1, 1) double {mustBeNonnegative} = 0.05
end

scores = double(scores(:));
scores = scores(isfinite(scores));
thres = [];
info = struct('x', [], 'y', [], 'extremaX', [], 'extremaY', [], 'bandwidth', []);
if numel(scores) < 3 || range(scores) == 0
    return;
end

% Silverman's rule, matching R's bw.nrd0, so results track scCancer's.
bw = opts.Bandwidth;
if bw == 0
    spread = min(std(scores), iqr(scores) / 1.349);
    if spread == 0, spread = std(scores); end
    if spread == 0, return; end
    bw = 0.9 * spread * numel(scores)^(-1/5);
end

% R's density() spans the data padded by three bandwidths; match it so the
% grid resolution and edge behaviour agree.
xi = linspace(min(scores) - 3*bw, max(scores) + 3*bw, opts.NumPoints);
y = ksdensity(scores, xi, 'Bandwidth', bw);
info.x = xi;
info.y = y;
info.bandwidth = bw;

% Local extrema are where the sign of the first difference turns over.
rising = diff(y) > 0;
extPos = find(diff(rising) ~= 0);
if isempty(extPos)
    return;
end
extY = y(extPos);
yMax = max(extY);

% A plateau shows up as a pair of extrema at nearly the same height. Drop
% both members of any such pair; they are ripples, not modes.
if numel(extPos) >= 3
    flat = abs(diff(extY)) < yMax * opts.MergeTol;
    drop = false(size(extY));
    drop([false, flat]) = true;
    drop([flat, false]) = true;
    extPos = extPos(~drop);
    extY = extY(~drop);
end
info.extremaX = xi(extPos);
info.extremaY = extY;
if numel(extPos) < 3
    return;
end

% Prominence of each extremum: how far it stands from its nearer neighbour,
% with the ends bounded by zero density.
padded = [0, extY, 0];
prominence = min(abs(padded(2:end-1) - padded(1:end-2)), ...
                 abs(padded(2:end-1) - padded(3:end)));

[~, byHeight] = sort(extY, 'descend');
tallest = byHeight(1);
second = byHeight(2);
if prominence(tallest) == 0 || ...
        prominence(second) / prominence(tallest) <= opts.ProminenceTol
    return;  % the second mode is a ripple on the shoulder of the first
end

% The valley is the lowest density between the two tallest peaks.
span = min(tallest, second):max(tallest, second);
[valleyY, k] = min(extY(span));

% Two peaks with no real trough between them are one peak with a ripple on
% it. Require the valley to sit a fair way below the shallower peak.
if 1 - valleyY / min(extY(tallest), extY(second)) < opts.MinDipDepth
    return;
end

thres = xi(extPos(span(k)));
end
