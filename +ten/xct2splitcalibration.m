function C = xct2splitcalibration(sce1, sce2, celltype1, celltype2, opts)
%XCT2SPLITCALIBRATION  Random-split null for TEN.XCT2NULLCONTROL's EmpP count.
%   C = TEN.XCT2SPLITCALIBRATION(SCE1, SCE2, CELLTYPE1, CELLTYPE2) pools
%   the two samples' CELLTYPE1 and CELLTYPE2 cells, splits each cell type
%   at random into two pseudo-samples of its original per-sample sizes, and
%   runs TEN.SCTENIFOLDXCT2 and TEN.XCT2NULLCONTROL on each split. Nothing
%   differs between the pseudo-samples, so the splits show how many pairs
%   reach EmpP <= 0.05 by chance with these cells, networks and settings.
%
%   WHY. TEN.XCT2NULLCONTROL's EmpP is not calibrated for a change between
%   conditions: real and shuffled networks are not exchangeable, and on six
%   random halves of one sample 17-43 pairs (mean 25) reached EmpP <= 0.05
%   against the 11.8 an exchangeable null predicts. A real comparison's
%   count is evidence of a difference only if it exceeds what its own random
%   splits give. See TEN.XCT2NULLCONTROL.
%
%   COST. Each split is one scTenifoldXct2 run plus NumShuffles control
%   runs: at the defaults, five splits take about NumSplits*21 times one
%   scTenifoldXct2 run (about an hour for 600 cells and 7,000 genes).
%   UseParallel runs the splits on the current parallel pool.
%
%   INPUTS:
%     sce1, sce2  - the two conditions, sharing one gene list
%     celltype1   - source cell type
%     celltype2   - target cell type
%
%   NAME-VALUE ARGUMENTS:
%     NumSplits   - random splits (default 5). The calibrated p cannot fall
%                   below 1/(NumSplits+1).
%     NumShuffles - control shuffles per split (default 20), as the real
%                   run's control used
%     XctArgs     - further TEN.SCTENIFOLDXCT2 arguments, as the real run's
%     RealCount   - the real comparison's NumEmpPBelow05 from
%                   TEN.XCT2NULLCONTROL; given, C says where it falls
%     Seed        - seed (default 1); split k uses Seed + k for both the
%                   split and its shuffles
%     UseParallel - run the splits with PARFOR (default false)
%     Verbose     - print the result (default true)
%
%   OUTPUT C, a struct:
%     NumSplits, NumShuffles, Seeds
%     SplitCounts  pairs at EmpP <= 0.05 in each split
%     SplitSigQ    pairs at q_value < 0.05 in each split (should be ~0)
%     MeanCount    mean of SplitCounts
%     RealCount    as given, or NaN
%     CalibratedP  (1 + #splits with a count >= RealCount)/(1 + NumSplits);
%                  NaN without RealCount
%     Summary      one sentence for the command line or a dialog
%
%   Example:
%       [T, grns] = ten.sctenifoldxct2(sce1, sce2, "Alpha cells", "Beta cells", false);
%       [Tn, S] = ten.xct2nullcontrol(sce1, sce2, "Alpha cells", "Beta cells", Grns=grns);
%       C = ten.xct2splitcalibration(sce1, sce2, "Alpha cells", "Beta cells", ...
%           RealCount=S.NumEmpPBelow05);
%
% see also: TEN.XCT2NULLCONTROL, TEN.SCTENIFOLDXCT2

arguments
    sce1 SingleCellExperiment
    sce2 SingleCellExperiment
    celltype1 (1,1) string
    celltype2 (1,1) string
    opts.NumSplits (1,1) double {mustBeInteger, mustBePositive} = 5
    opts.NumShuffles (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.XctArgs cell = {}
    opts.RealCount (1,1) double = NaN
    opts.Seed (1,1) double = 1
    opts.UseParallel (1,1) logical = false
    opts.Verbose (1,1) logical = true
end

if ~isequal(string(sce1.g(:)), string(sce2.g(:)))
    error("TEN:XCT2SPLITCALIBRATION:GeneMismatch", ...
        "SCE1 and SCE2 must share one gene list, as TEN.SCTENIFOLDXCT2 requires.");
end
types = [celltype1, celltype2];
in1 = ismember(string(sce1.c_cell_type_tx), types);
in2 = ismember(string(sce2.c_cell_type_tx), types);
X = [sce1.X(:, in1), sce2.X(:, in2)];
ct = [string(sce1.c_cell_type_tx(in1)); string(sce2.c_cell_type_tx(in2))];
g = sce1.g;
% Per cell type, how many of its cells sample 1 had: each split keeps it.
num1 = [nnz(string(sce1.c_cell_type_tx(in1)) == celltype1), ...
    nnz(string(sce1.c_cell_type_tx(in1)) == celltype2)];

n = opts.NumSplits;
seeds = opts.Seed + (1:n);
counts = nan(1, n);
sigq = nan(1, n);
xctArgs = opts.XctArgs;
numShuffles = opts.NumShuffles;
if opts.UseParallel
    parfor k = 1:n
        [counts(k), sigq(k)] = i_onesplit(X, g, ct, types, num1, seeds(k), ...
            xctArgs, numShuffles);
    end
else
    for k = 1:n
        if opts.Verbose
            fprintf('[xct2split] split %d of %d...\n', k, n);
        end
        [counts(k), sigq(k)] = i_onesplit(X, g, ct, types, num1, seeds(k), ...
            xctArgs, numShuffles);
    end
end

realCount = opts.RealCount;
calibratedP = NaN;
summary = sprintf(['%d random splits of the same cells gave %s pairs at ' ...
    'EmpP <= 0.05 (mean %.1f)'], n, strjoin(string(counts), ", "), mean(counts));
if ~isnan(realCount)
    calibratedP = (1 + nnz(counts >= realCount))/(1 + n);
    summary = sprintf(['%s; the real comparison''s %d is exceeded or matched ' ...
        'by %d of them (calibrated p = %.3g).'], summary, realCount, ...
        nnz(counts >= realCount), calibratedP);
else
    summary = [summary, '.'];
end
if opts.Verbose
    fprintf('[xct2split] %s\n', summary);
end

C = struct('NumSplits', n, 'NumShuffles', numShuffles, 'Seeds', seeds, ...
    'SplitCounts', counts, 'SplitSigQ', sigq, 'MeanCount', mean(counts), ...
    'RealCount', realCount, 'CalibratedP', calibratedP, ...
    'Summary', string(summary));
end


function [count, sigq] = i_onesplit(X, g, ct, types, num1, seed, xctArgs, numShuffles) %#ok<INUSD> xctArgs and numShuffles are used inside evalc
% One random split: each cell type's pooled cells reassigned to two
% pseudo-samples of the original sizes, then scTenifoldXct2 and its control.
rs = RandStream("twister", "Seed", seed);
inA = false(numel(ct), 1);
for t = 1:2
    idx = find(ct == types(t));
    pick = randperm(rs, numel(idx), num1(t));
    inA(idx(pick)) = true;
end
sceA = SingleCellExperiment(X(:, inA), g);
sceA.c_cell_type_tx = ct(inA);
sceB = SingleCellExperiment(X(:, ~inA), g);
sceB.c_cell_type_tx = ct(~inA);
[~, T, grns] = evalc("ten.sctenifoldxct2(sceA, sceB, types(1), types(2), false, xctArgs{:})"); %#ok<ASGLU> grns is used inside the next evalc
sigq = nnz(T.q_value < 0.05);
[~, ~, S] = evalc("ten.xct2nullcontrol(sceA, sceB, types(1), types(2), Grns=grns, XctArgs=xctArgs, NumShuffles=numShuffles, Seed=seed, Verbose=false)");
count = S.NumEmpPBelow05;
end
