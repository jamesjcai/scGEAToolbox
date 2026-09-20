function [isDoublet, score, info] = sc_scrublet(X, opts)
%SC_SCRUBLET Detect doublets by comparison with simulated doublets.
%
%   [isDoublet, score] = SC_SCRUBLET(X) scores every cell by how much of its
%   neighbourhood is made of artificial doublets, and calls the cells above
%   an automatic threshold.
%
%   This is a native MATLAB reimplementation of Scrublet (Wolock, Lopez &
%   Klein, Cell Systems 2019, PMID:30954476). No Python is involved.
%
%   The idea is to answer "what would a doublet look like in this dataset?"
%   with the dataset itself. Pairs of observed cells are added together to
%   make synthetic doublets, the real and synthetic cells are embedded
%   together, and each real cell is scored by the fraction of its nearest
%   neighbours that are synthetic. A real cell sitting in a cloud of
%   simulated doublets probably is one.
%
%   The method only sees doublets formed from two *different* cell types.
%   A doublet of two cells of the same type lands on top of that type and
%   is invisible here, which is why the simulated scores are bimodal and
%   why the detected rate is always below the true rate.
%
%   USAGE:
%     [isDoublet, score] = sc_scrublet(sce.X);
%     sce.setCellAttribute('doublet_score', score);
%     sce = sce.selectcells(~isDoublet);
%
%   INPUT:
%     X - nGenes-by-nCells matrix of raw UMI counts, full or sparse. Raw
%         counts, not normalised: the simulation adds transcriptomes, and
%         that is only meaningful on counts.
%
%   NAME-VALUE:
%     SimDoubletRatio     - synthetic doublets to simulate, as a multiple of
%                           the number of cells (default 2).
%     NumNeighbors        - neighbourhood size k (default
%                           round(0.5*sqrt(nCells))). The search itself uses
%                           k*(1 + nSim/nCells) neighbours, so that k of
%                           them are expected to be real cells.
%     ExpectedDoubletRate - prior doublet rate of the experiment
%                           (default 0.1). Roughly 0.008 per 1000 cells
%                           loaded on a 10x run.
%     StdevDoubletRate    - uncertainty on that prior (default 0.02), used
%                           only for the standard error in INFO.
%     NumPCs              - principal components of the joint embedding
%                           (default 30).
%     MinCounts, MinCells - a gene is kept when at least MinCells cells have
%                           at least MinCounts normalised counts of it
%                           (defaults 3 and 3).
%     HVGPercentile       - percentile of the variability statistic above
%                           which genes are kept (default 85, i.e. the top
%                           15%). Set [] to keep every gene that passes the
%                           count filter.
%     Threshold           - call threshold on the score. Default [] picks
%                           the valley between the two modes of the
%                           simulated-doublet scores.
%     Seed                - RNG seed (default [], no reseed).
%     Verbose             - print progress (default false).
%
%   OUTPUTS:
%     isDoublet - nCells-by-1 logical, true for called doublets.
%     score     - nCells-by-1 doublet score in [0, 1].
%     info      - struct with the simulated-doublet scores (simScore), the
%                 threshold used and whether it was found automatically,
%                 the parent cell indices of each simulated doublet, the
%                 neighbourhood sizes, the standard error of each score,
%                 and the detected and estimated doublet rates.
%                 detectableFraction is the share of simulated doublets that
%                 clear the threshold, i.e. how much power the method has on
%                 this dataset; well-separated data sits near 0.75, and a low
%                 value raises a warning.
%
%   NOTE: the gene filter uses this toolbox's own variability statistic
%   (SC_HVG, CV^2 against a fitted mean-variance trend) in place of
%   Scrublet's v-score. Both measure above-Poisson variability and pick out
%   the same kind of gene; the scores are not numerically identical to the
%   Python package's.
%
%   See also SC_HVG, SC_QCFILTER, GUI.CALLBACK_DOUBLETDETECTION.

arguments
    X {mustBeNumeric, mustBeNonempty}
    opts.SimDoubletRatio (1,1) double {mustBePositive} = 2
    opts.NumNeighbors double {mustBeScalarOrEmpty} = []
    opts.ExpectedDoubletRate (1,1) double {mustBePositive} = 0.1
    opts.StdevDoubletRate (1,1) double {mustBeNonnegative} = 0.02
    opts.NumPCs (1,1) double {mustBePositive, mustBeInteger} = 30
    opts.MinCounts (1,1) double {mustBeNonnegative} = 3
    opts.MinCells (1,1) double {mustBeNonnegative} = 3
    opts.HVGPercentile double {mustBeScalarOrEmpty} = 85
    opts.Threshold double {mustBeScalarOrEmpty} = []
    opts.Seed = []
    opts.Verbose (1,1) logical = false
end

[numGenes, numCells] = size(X);
if numCells < 10
    error("sc_scrublet:tooFewCells", ...
        "Doublet detection needs more than %d cells.", numCells);
end
if opts.ExpectedDoubletRate >= 1
    error("sc_scrublet:rateOutOfRange", ...
        "EXPECTEDDOUBLETRATE is a proportion; %g is not one.", ...
        opts.ExpectedDoubletRate);
end

if ~isempty(opts.Seed)
    rng(opts.Seed);
end

% ---- Simulate doublets by adding pairs of observed cells ---------------
numSim = max(1, round(opts.SimDoubletRatio*numCells));
parents = randi(numCells, numSim, 2);
Xsim = X(:, parents(:, 1)) + X(:, parents(:, 2));
if opts.Verbose
    fprintf("sc_scrublet: %d cells, %d simulated doublets\n", numCells, numSim);
end

% ---- Normalise both to the same library size --------------------------
% The target is the mean observed library size, so the simulated doublets
% end up on the scale of real cells rather than at twice it.
totalObs = full(sum(X, 1));
totalSim = full(sum(Xsim, 1));
target = mean(totalObs);
Xobs = i_libnorm(X, totalObs, target);
Xsim = i_libnorm(Xsim, totalSim, target);

% ---- Gene filter, decided on the observed cells only -------------------
% Nothing about the simulation should be allowed to choose the genes: the
% synthetic cells are meant to be projected into the real cells' space.
keep = sum(Xobs >= opts.MinCounts, 2) >= opts.MinCells;
if ~any(keep)
    error("sc_scrublet:noGenes", ...
        ['No gene has at least %g counts in at least %d cells. ' ...
        'Lower MINCOUNTS or MINCELLS, or check that X holds raw counts.'], ...
        opts.MinCounts, opts.MinCells);
end
if ~isempty(opts.HVGPercentile) && sum(keep) > 50
    idxKeep = find(keep);
    T = sc_hvg(Xobs(idxKeep, :), [], false, false, false);
    cutoff = prctile(T.fitratio, opts.HVGPercentile);
    isVariable = T.fitratio >= cutoff & isfinite(T.fitratio);
    if sum(isVariable) >= 10
        keep = false(numGenes, 1);
        keep(idxKeep(isVariable)) = true;
    end
end
Xobs = Xobs(keep, :);
Xsim = Xsim(keep, :);
if opts.Verbose
    fprintf("sc_scrublet: %d of %d genes retained\n", sum(keep), numGenes);
end

% ---- Standardise, with the observed cells' mean and sd ------------------
geneMean = full(mean(Xobs, 2));
geneStd = full(std(Xobs, 0, 2));
geneStd(geneStd < eps) = 1;
Xobs = single(full((Xobs - geneMean)./geneStd)).';
Xsim = single(full((Xsim - geneMean)./geneStd)).';

% ---- PCA fitted on the observed cells, simulated cells projected -------
numPCs = min([opts.NumPCs, size(Xobs, 2), numCells - 1]);
coeff = pca(Xobs, "NumComponents", numPCs, "Centered", false, ...
    "Algorithm", "eig");
manifold = double([Xobs*coeff; Xsim*coeff]);
clear Xobs Xsim;

% ---- Nearest neighbours over the real and simulated cells together -----
k = opts.NumNeighbors;
if isempty(k)
    k = round(0.5*sqrt(numCells));
end
k = max(1, round(k));
ratio = numSim/numCells;
kAdj = min(round(k*(1 + ratio)), size(manifold, 1) - 1);

nbr = knnsearch(manifold, manifold, "K", kAdj + 1);
nbr = nbr(:, 2:end);                    % drop each point's own index
numSimNeighbors = sum(nbr > numCells, 2);

% ---- Score: posterior probability that a cell is a doublet -------------
% From Bayes' rule on "is this neighbourhood a doublet neighbourhood",
% with the simulated cells standing in for the doublet class and the
% experiment's prior rate correcting for how many were simulated.
rho = opts.ExpectedDoubletRate;
q = (numSimNeighbors + 1)/(kAdj + 2);
denom = 1 - rho - q.*(1 - rho - rho/ratio);
allScore = q*rho/ratio./denom;

seQ = sqrt(q.*(1 - q)/(kAdj + 2));
seRho = opts.StdevDoubletRate;
allSe = q*rho/ratio./denom.^2 .* ...
    sqrt((seQ./q*(1 - rho)).^2 + (seRho/rho*(1 - q)).^2);

score = allScore(1:numCells);
simScore = allScore(numCells+1:end);

% ---- Threshold ---------------------------------------------------------
% The simulated scores are bimodal: doublets of two different types score
% high, doublets of one type land on that type and score low. The valley
% between the modes separates the detectable doublets from the rest.
threshold = opts.Threshold;
isAuto = isempty(threshold);
if isAuto
    threshold = i_minimumthreshold(simScore);
end
if isempty(threshold)
    warning("sc_scrublet:noThreshold", ...
        ['The simulated doublet scores are not bimodal, so no threshold ' ...
        'could be found. Inspect INFO.simScore and pass THRESHOLD.']);
    isDoublet = false(numCells, 1);
    threshold = NaN;
else
    isDoublet = score > threshold;
end

detectedRate = mean(isDoublet);
detectableFraction = mean(simScore > threshold);

% How many of the simulated doublets clear the threshold says how much power
% the method has on this dataset. When most of them look like singlets there
% is no real valley in the histogram, the automatic threshold is cutting
% noise, and the calls should not be trusted. That is what one homogeneous
% population looks like; datasets with distinct types sit near 0.75.
if isAuto && isfinite(threshold) && detectableFraction < 0.3
    warning("sc_scrublet:weakSeparation", ...
        ['Only %.0f%% of simulated doublets score above the automatic ' ...
        'threshold, so most of them are indistinguishable from single ' ...
        'cells and the threshold rests on a weak split. Inspect ' ...
        'INFO.simScore and set THRESHOLD explicitly before trusting ' ...
        'these calls.'], 100*detectableFraction);
end
info = struct( ...
    "simScore", simScore, ...
    "threshold", threshold, ...
    "thresholdIsAutomatic", isAuto, ...
    "doubletParents", parents, ...
    "numNeighbors", k, ...
    "numNeighborsAdjusted", kAdj, ...
    "numSimulated", numSim, ...
    "standardError", allSe(1:numCells), ...
    "detectedDoubletRate", detectedRate, ...
    "detectableFraction", detectableFraction, ...
    "overallDoubletRate", detectedRate/max(detectableFraction, eps), ...
    "genesUsed", sum(keep));

if opts.Verbose
    fprintf("sc_scrublet: threshold %.4f, %s called (%.1f%%)\n", ...
        threshold, pkg.i_plural(sum(isDoublet), 'doublet'), 100*detectedRate);
    fprintf("sc_scrublet: %.1f%% of simulated doublets are detectable, " + ...
        "so the true rate is nearer %.1f%%\n", ...
        100*detectableFraction, 100*info.overallDoubletRate);
end

end


function Xn = i_libnorm(Xn, total, target)
% Scale each cell to TARGET total counts. Empty cells are left at zero
% rather than turned into NaN.
total(total == 0) = 1;
Xn = Xn.*(target./total);
end


function threshold = i_minimumthreshold(v)
% The "minimum" threshold of Prewitt & Mendelsohn, as used by Scrublet
% through skimage: smooth the histogram until exactly two maxima remain,
% then take the lowest point between them. Returns [] when the histogram
% never settles at two modes, which means the scores are not bimodal and
% there is no valley to find.

threshold = [];
v = v(isfinite(v));
if numel(v) < 10 || range(v) < eps
    return;
end

numBins = 256;
edges = linspace(min(v), max(v), numBins + 1);
centres = edges(1:end-1) + diff(edges)/2;
counts = double(histcounts(v, edges));

peaks = [];
for iter = 1:10000
    counts = i_smooth3(counts);
    peaks = i_localmaxima(counts);
    if numel(peaks) < 3
        break;
    end
end
if numel(peaks) ~= 2
    return;
end

[~, offset] = min(counts(peaks(1):peaks(2)));
threshold = centres(peaks(1) + offset - 1);
end


function counts = i_smooth3(counts)
% Three-point running mean with the edges reflected, matching the uniform
% filter skimage applies. MOVMEAN shrinks its window at the ends instead,
% which changes whether the first and last bins come out as maxima.
padded = [counts(1), counts, counts(end)];
counts = (padded(1:end-2) + padded(2:end-1) + padded(3:end))/3;
end


function idx = i_localmaxima(counts)
% A rising-edge detector, matching skimage's threshold_minimum: a bin is a
% maximum when the histogram stops rising at it. Unlike ISLOCALMAX this
% counts the first bin when the histogram starts by falling -- which is the
% usual shape here, since the low-scoring mode of the simulated doublets
% piles up against zero. Treating that mode as no maximum at all leaves only
% interior bumps to pair up, and the valley lands far to the right of the
% real one, above most of the doublets it was supposed to separate.
idx = zeros(1, numel(counts));
found = 0;
rising = true;
for i = 1:numel(counts)-1
    if rising
        if counts(i+1) < counts(i)
            rising = false;
            found = found + 1;
            idx(found) = i;
        end
    elseif counts(i+1) > counts(i)
        rising = true;
    end
end
idx = idx(1:found);
end
