function [stats, M] = ml_mine(x, y, opts)
%ML_MINE Maximal information-based nonparametric exploration statistics.
%   STATS = run.ml_mine(X, Y) returns MIC, MAS, MEV, MCN, MCN_GENERAL and
%   TIC for the paired samples X and Y, computed with the ApproxMaxMI
%   heuristic of Reshef et al. (Science 2011, Supplementary Online
%   Material, Algorithms 1-5).
%
%   [STATS, M] = run.ml_mine(X, Y) also returns the characteristic
%   matrix. M(i, j) is the normalised mutual information of the best grid
%   with i + 1 bins on Y and j + 1 bins on X, and is zero where the grid
%   exceeds the (i + 1)(j + 1) <= B limit.
%
%   Name-value arguments:
%     Alpha - grid-size exponent: B = max(n^Alpha, 4) for Alpha in (0, 1],
%             or B = Alpha directly when Alpha >= 4. Default 0.6.
%     C     - clump factor: at most C*x clumps are kept when optimising
%             x columns. Default 15.
%
%   This is an independent implementation written from the published
%   pseudocode, not a translation of minepy (GPL-3). It reproduces
%   minepy's est="mic_approx" output -- statistics and characteristic
%   matrix -- to about 1e-13, including on heavily tied data; see
%   tests/mlMineTest.m for the three conventions that required. Speed is
%   close to minepy's C: about 0.13 s a pair at 2000 samples.
%
%   See also run.mex_minepy.

arguments
    x (:, 1) double {mustBeReal, mustBeFinite}
    y (:, 1) double {mustBeReal, mustBeFinite}
    opts.Alpha (1, 1) double {mustBePositive} = 0.6
    opts.C (1, 1) double {mustBePositive} = 15
end

n = numel(x);
if numel(y) ~= n
    error("ml_mine:sizeMismatch", "X and Y must have the same number of elements.");
end
% minepy returns MIC = 1 from two points and fails outright on one, so
% below four samples there is nothing coherent to match.
if n < 4
    error("ml_mine:tooFewSamples", "MINE needs at least 4 samples; got %d.", n);
end
if opts.Alpha > 1 && opts.Alpha < 4
    error("ml_mine:badAlpha", "Alpha must be in (0, 1] or >= 4.");
end

if opts.Alpha <= 1
    B = max(n^opts.Alpha, 4);
else
    B = opts.Alpha;
end
B = min(B, n);
maxBins = floor(B/2);

% Ixy(a, b): best normalised MI with a columns on X and b rows on Y, from
% each orientation of the search. The two orientations are combined
% afterwards.
Ixy = zeros(maxBins);
Iyx = zeros(maxBins);
xAxis = sortAxis(x);
yAxis = sortAxis(y);
for rowBins = 2:maxBins
    colBins = floor(B/rowBins);
    Ixy(2:colBins, rowBins) = approxMaxMI(xAxis, yAxis, colBins, rowBins, opts.C*colBins);
    Iyx(rowBins, 2:colBins) = approxMaxMI(yAxis, xAxis, colBins, rowBins, opts.C*colBins);
end

% Characteristic matrix indexed (Y bins, X bins).
[yBins, xBins] = ndgrid(2:maxBins);
isValid = (xBins.*yBins) <= B;
Ibest = max(Ixy(2:end, 2:end), Iyx(2:end, 2:end)).';
M = zeros(maxBins - 1);
M(isValid) = Ibest(isValid);

stats = summariseMatrix(M, xBins, yBins, isValid);
end

function stats = summariseMatrix(M, xBins, yBins, isValid)
mic = max(M(isValid));

Mt = M.';
bothValid = isValid & isValid.';
mas = max(abs(M(bothValid) - Mt(bothValid)), [], "all");
if isempty(mas)
    mas = 0;
end

onEdge = isValid & ((xBins == 2) | (yBins == 2));
mev = max(M(onEdge));

% MCN counts a grid as reaching MIC when it comes within 1e-4 of it. That
% is minepy's convention: fitting it to minepy's output bounds the absolute
% tolerance to [9.1e-5, 1.08e-4], and no relative tolerance fits at all.
tolerance = 1e-4;
gridSize = log2(xBins.*yBins);
mcn = min(gridSize(isValid & (M >= mic - tolerance)));
mcnGeneral = min(gridSize(isValid & (M >= mic*mic - tolerance)));

stats = struct("mic", mic, "mas", mas, "mev", mev, "mcn", mcn, ...
    "mcn_general", mcnGeneral, "tic", sum(M(isValid)));
end

function sortedAxis = sortAxis(v)
% Sort order and tie runs of one variable. They do not depend on the grid,
% so they are computed once rather than on every ApproxMaxMI call.
[vs, order] = sort(v);
isRunStart = [true; diff(vs) ~= 0];
runStart = find(isRunStart);
sortedAxis = struct("order", order, "runId", cumsum(isRunStart), ...
    "runStart", runStart, "runEnd", [runStart(2:end) - 1; numel(v)]);
end

function I = approxMaxMI(xAxis, yAxis, colBins, rowBins, maxClumps)
% Best normalised MI for 2..colBins columns on the x axis, with the y axis
% equipartitioned into ROWBINS rows (SOM Algorithm 2).
Q = equipartitionAxis(yAxis, rowBins);
Qs = Q(xAxis.order);
P = superclumpsPartition(xAxis.runId, Qs, maxClumps);
I = optimizeXAxis(Qs, P, colBins);
end

function Q = equipartitionAxis(sortedAxis, numBins)
% Assign each point to one of NUMBINS rows of near-equal size, keeping
% tied values together (SOM Algorithm 3, EquipartitionYAxis).
% The greedy pass is sequential, so it stays a loop; it assigns a row to
% each tie run and expands to points afterwards.
runStart = sortedAxis.runStart;
runSize = sortedAxis.runEnd - runStart + 1;
numRuns = numel(runStart);
n = numel(sortedAxis.order);
runRow = zeros(numRuns, 1);
currentRow = 1;
rowCount = 0;
desiredSize = n/numBins;
for r = 1:numRuns
    s = runSize(r);
    if (rowCount ~= 0) && (abs(rowCount + s - desiredSize) >= abs(rowCount - desiredSize))
        currentRow = currentRow + 1;
        rowCount = 0;
        desiredSize = (n - runStart(r) + 1)/(numBins - currentRow + 1);
    end
    runRow(r) = currentRow;
    rowCount = rowCount + s;
end
Q = zeros(n, 1);
Q(sortedAxis.order) = runRow(sortedAxis.runId);
end

function P = superclumpsPartition(runId, Qs, maxClumps)
% Clump label for each point, points sorted by x (SOM Algorithm 4). Too
% many clumps are merged by equipartitioning the clump labels themselves,
% which keeps every clump whole.
P = clumpsPartition(runId, Qs);
if P(end) > maxClumps
    P = equipartitionAxis(sortAxis(P), floor(maxClumps));
end
end

function P = clumpsPartition(runId, Qs)
% Consecutive points in the same row form a clump. Points sharing an x
% value (one RUNID) cannot be separated by a column, so a tie spanning
% several rows becomes a clump of its own.
isMixed = accumarray(runId, Qs, [], @max) ~= accumarray(runId, Qs, [], @min);
label = Qs;
mixedPoint = isMixed(runId);
label(mixedPoint) = -runId(mixedPoint);
P = cumsum([true; diff(label) ~= 0]);
end

function I = optimizeXAxis(Qs, P, colBins)
% Exact dynamic programme over clump boundaries (SOM Algorithm 5). The
% MI of a column partition is H(Q) - sum_j (n_j/n) H(Q | column j), and
% the sum is additive over columns, so G(l, t) = the least
% sum_j n_j H(Q | column j) over l columns covering clumps 1..t.
n = numel(Qs);
numClumps = P(end);
numRows = max(Qs);

counts = accumarray([P, Qs], 1, [numClumps, numRows]);
cumCounts = [zeros(1, numRows); cumsum(counts, 1)];
cumTotal = sum(cumCounts, 2);

% cost(s + 1, t + 1): n_col * H(Q | col) for a column spanning clumps s+1..t.
cost = xlogx(cumTotal.' - cumTotal);
for r = 1:numRows
    cost = cost - xlogx(cumCounts(:, r).' - cumCounts(:, r));
end
cost(tril(true(numClumps + 1))) = Inf;

rowProb = counts;
rowProb = sum(rowProb, 1)/n;
HQ = -sum(rowProb(rowProb > 0).*log(rowProb(rowProb > 0)));

% Normalise by log(min(columns, rows)) with the rows the equipartition
% actually produced, which ties can leave fewer than asked for. The
% column count stays nominal even when there are fewer clumps. Both
% choices are minepy's, established by comparing against its output.
I = zeros(colBins - 1, 1);
G = cost(1, :);
best = -Inf;
for l = 2:colBins
    % Past one column per clump no partition can improve, so the DP
    % step is skipped and BEST carries forward.
    % l columns need at least l clumps, so only boundaries s >= l - 1
    % and t >= l are live (indices shifted by one for s = 0).
    if l <= numClumps
        live = l:numClumps;
        G(live + 1) = min(G(live).' + cost(live, live + 1), [], 1);
        G(1:l) = Inf;
        best = max(best, HQ - G(end)/n);
    end
    usableBins = min(l, numRows);
    if usableBins > 1
        I(l - 1) = max(best, 0)/log(usableBins);
    end
end
end

function v = xlogx(m)
% m*log(m) with 0*log(0) = 0. M holds counts -- non-negative integers --
% so clamping at 1 changes nothing but the zeros, and avoids masking.
v = m.*log(max(m, 1));
end
