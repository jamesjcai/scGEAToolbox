function [segmented, breakpoints, scaleSeq] = e_cnvsegment(Y, chr, opts)
%E_CNVSEGMENT  Piecewise-constant segmentation of a CNV profile.
%
%   [segmented, breakpoints] = PKG.E_CNVSEGMENT(Y, chr) turns a per-gene
%   copy-number profile into piecewise-constant segments by merging adjacent
%   regions cheapest-first until the cheapest merge left would cross a real
%   boundary. This is the step copykat and SCEVAN add on top of the
%   inferCNV-style moving average: it localises breakpoints and replaces a
%   smooth ridge with the discrete gains and losses a karyotype is made of.
%
%   The segmentation is joint across cells: one set of breakpoints is fitted
%   to all of them at once. Cells of a clone share their boundaries, so
%   pooling finds those boundaries far more sharply than segmenting each
%   cell alone, and every cell ends up on the same coordinates, which is
%   what makes the columns comparable afterwards.
%
%   INPUTS:
%     Y   - nGenes-by-nCells profile, rows ordered along the genome.
%     chr - nGenes-by-1 chromosome index per row. Merges never cross a
%           chromosome, so each one is segmented independently.
%
%   NAME-VALUE:
%     StopRatio - stop merging once the cheapest merge costs more than this
%                 multiple of the noise level (default 2). The cost is
%                 normalised so that ~1 is what an uninteresting boundary
%                 costs, which is what makes this number transferable.
%     MaxMerges - safety cap on the number of merges (default Inf).
%     UseSingle - work in single precision (default true).
%
%   OUTPUTS:
%     segmented   - nGenes-by-nCells, piecewise constant within a segment.
%     breakpoints - indices where a new segment starts, plus a final
%                   sentinel one past the end.
%     scaleSeq    - the merge cost at each step, in order.
%
%   Extracted from RUN.ML_SCEVAN so that SC_INFERCNV can offer the same
%   segmentation without a second copy of it. See also SC_INFERCNV,
%   SC_MALIGNSCORE, RUN.ML_SCEVAN.

arguments
    Y {mustBeNumeric, mustBeNonempty}
    chr (:, 1) {mustBeNumeric}
    opts.StopRatio (1, 1) double {mustBePositive} = 2
    opts.MaxMerges (1, 1) double {mustBePositive} = inf
    opts.UseSingle (1, 1) logical = true
end

if numel(chr) ~= size(Y, 1)
    error("e_cnvsegment:chrSize", ...
        "CHR has %d entries but Y has %d rows.", numel(chr), size(Y, 1));
end
opts.segStopRatio = opts.StopRatio;
opts.maxMerges = opts.MaxMerges;
opts.segUseSingle = opts.UseSingle;
[segmented, breakpoints, scaleSeq] = i_segment(Y, chr, opts);
end

% =========================================================================
function [Useg, breakpoints, scaleSeq] = i_segment(Y, chr, opts)
% Greedy piecewise-constant segmentation along the genome (the
% Mumford-Shah special case of [1]). Adjacent regions are merged
% cheapest-first until the cheapest merge left would cross a real boundary.
%
% The cost of merging regions i and j is
%
%     cost = (li*lj/(li+lj)) * sum_over_cells (ui - uj)^2                [1]
%
% Under the hypothesis that the two regions share one true value, the
% difference of their means has variance sigma^2*(1/li + 1/lj) per cell,
% and li*lj/(li+lj) is precisely the reciprocal of that factor. So
%
%     ratio = cost / (nCells * sigma^2)
%
% is a chi-square statistic on nCells degrees of freedom divided by its own
% degrees of freedom: about 1 when the boundary is noise, and growing
% without limit at a real one. Merging stops once the cheapest available
% merge exceeds opts.segStopRatio.
%
% Normalising this way, rather than reducing the cost to a plain
% root-mean-square difference, is what makes the rule work. The region
% sizes must stay in: two single genes that differ are unremarkable, while
% two 500-gene regions differing by the same amount are a breakpoint, and
% only the li*lj/(li+lj) factor distinguishes them.
%
% THREE FIXES over the original implementation:
%
%   1. The stopping test was inverted. It stopped when the increase in
%      lambda was SMALL -- which is exactly what happens at the first few
%      merges -- so segmentation halted after a handful of them and left
%      CNA at gene resolution. It also compared a lambda increment, in cost
%      units, against nu, in expression units; the two are not commensurate
%      at all, which is why no value of the old betaStop gave sensible
%      behaviour. opts.betaStop is therefore replaced by
%      opts.segStopRatio, whose scale is meaningful.
%   2. Regions could merge across chromosome boundaries, smearing an event
%      on one chromosome into the start of the next. Those merges are now
%      forbidden outright.
%   3. Regions were stored as rows of a matrix and deleted one per merge.
%      That is O(nRegions^2 * nCells) of pure copying and only went
%      unnoticed because bug 1 stopped the loop almost immediately; with it
%      fixed, a 2800-region by 8000-cell problem spent ~40 s shuffling
%      memory. Regions now live in fixed slots joined by a linked list, so
%      a merge costs O(nCells).

[nGenes, nCells] = size(Y);
if opts.segUseSingle && ~isa(Y, 'single')
    Y = single(Y);
end
chr = double(chr(:));

% Per-gene noise SD, from adjacent differences. For d ~ N(0, 2*sigma^2),
% median|d| = 0.6745*sqrt(2)*sigma.
sigma = total_variability(Y) / (0.6745 * sqrt(2));
denom = nCells * sigma^2;

Umean = Y;
len = ones(nGenes, 1);
nxt = [2:nGenes, 0].';
prv = [0, 1:nGenes-1].';
active = true(nGenes, 1);

% Initial costs: every region is one gene, so effLen = 1*1/(1+1) = 1/2.
cost = inf(nGenes, 1);
cost(1:nGenes-1) = 0.5 * sum(double(Y(1:end-1, :) - Y(2:end, :)).^2, 2);
cost([chr(1:end-1) ~= chr(2:end); false]) = Inf;

maxMerges = min(opts.maxMerges, nGenes - 1);
scaleSeq = zeros(maxMerges, 1);
mergesDone = 0;

while mergesDone < maxMerges
    [minCost, i] = min(cost);
    if ~isfinite(minCost)
        break;   % every remaining boundary is a chromosome boundary
    end
    j = nxt(i);

    ratio = double(minCost) / denom;
    if ratio > opts.segStopRatio
        break;   % the cheapest merge left joins genuinely different regions
    end

    newLen = len(i) + len(j);
    Umean(i, :) = (len(i)*Umean(i, :) + len(j)*Umean(j, :)) / newLen;
    len(i) = newLen;

    % Unlink j and repair the costs of the two boundaries it touched.
    active(j) = false;
    cost(j) = Inf;
    nxt(i) = nxt(j);
    if nxt(i) > 0
        prv(nxt(i)) = i;
    end
    cost(i) = pair_cost(Umean, len, chr, i, nxt(i));
    if prv(i) > 0
        cost(prv(i)) = pair_cost(Umean, len, chr, prv(i), i);
    end

    mergesDone = mergesDone + 1;
    scaleSeq(mergesDone) = ratio;
end
scaleSeq = scaleSeq(1:mergesDone);

[Useg, breakpoints] = expand_regions(Umean, len, active, nxt);
end

function c = pair_cost(Umean, len, chr, i, j)
% Cost of merging region i with the region that follows it. Infinite when
% there is no next region or when it sits on another chromosome.
if j <= 0 || chr(i) ~= chr(j)
    c = Inf;
    return;
end
dv = double(Umean(i, :) - Umean(j, :));
c = (len(i)*len(j) / (len(i) + len(j))) * sum(dv.^2);
end

function nu = total_variability(Y)
% Simple proxy for "total variability" nu [1].
% Use median absolute adjacent difference across genome, aggregated across cells.
%
% Ties are excluded. Any profile that has been flattened somewhere -- and
% SC_INFERCNV flattens the whole dead zone its reference cells occupy, which
% on real input is most of the genome -- carries a large block of exactly
% equal neighbours. Those are an artefact of the flattening and say nothing
% about the noise, but they are numerous enough to take the median to zero,
% and a zero noise level makes every candidate merge look infinitely
% expensive: segmentation then stops before it starts and returns one
% segment per gene.
d = abs(double(reshape(diff(Y, 1, 1), [], 1)));
d = d(d > 0);
if isempty(d)
    nu = 1e-6;
    return;
end
nu = median(d);
if ~isfinite(nu) || nu == 0
    nu = 1e-6;
end
end

function [Useg, breakpoints] = expand_regions(Umean, len, active, nxt)
% Expand region means back to a gene-level piecewise-constant signal by
% walking the linked list of surviving regions.
nR = sum(active);
nCells = size(Umean, 2);
nGenes = sum(len(active));

Useg = zeros(nGenes, nCells, 'like', Umean);
breakpoints = zeros(nR + 1, 1);
breakpoints(1) = 1;

idx = find(active, 1);
p = 1;
r = 0;
while idx > 0
    r = r + 1;
    L = len(idx);
    Useg(p:p+L-1, :) = repmat(Umean(idx, :), L, 1);
    p = p + L;
    breakpoints(r+1) = p;
    idx = nxt(idx);
end
% breakpoints are 1..nGenes+1 style; last equals nGenes+1
end
