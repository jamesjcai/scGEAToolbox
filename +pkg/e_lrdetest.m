function [idx, pval] = e_lrdetest(X, y, k, progressFcn)
%E_LRDETEST  Rank genes by a likelihood-ratio test of group membership.
%
%   idx = PKG.E_LRDETEST(X, y, k) returns the indices of the K genes most
%   strongly associated with the two-group label Y, judged by a likelihood
%   ratio test of a logistic regression on the gene against an
%   intercept-only model.
%
%   [idx, pval] = PKG.E_LRDETEST(...) also returns the per-gene p-values,
%   with Inf for genes that were not testable.
%
%   PKG.E_LRDETEST(X, y, k, progressFcn) calls progressFcn(fraction) as it
%   goes, for a caller that wants to drive a waitbar.
%
%   X is CELLS-by-GENES, matching the orientation the caller gets from
%   sce.X'. Y is one label per cell, 0 or 1.
%
%   This was a nested function inside +gui/callback_Brush4Markers, where it
%   was written as
%
%       n = size(X, 1);              % number of CELLS
%       p_val = zeros(n, 1);         % but holds per-GENE p-values
%       for x = 1:size(X, 1)         % iterates over CELLS
%           model_data = table(X(:, x), ...);   % indexes a GENE
%
%   so it tested genes 1..nCells and never looked at the rest, while the
%   caller used the returned indices to pick gene names. Measured on 40
%   cells and 120 genes with five planted markers at genes 100-104: not one
%   of the five was ever tested and the callback reported five arbitrary
%   noise genes as the markers of the brushed selection. With fewer genes
%   than cells it raised MATLAB:badsubscript instead. Moving it here makes
%   the arithmetic testable without a figure.
%
%   See also GUI.CALLBACK_BRUSH4MARKERS, SC_DEG.

arguments
    X {mustBeNumeric, mustBeReal}
    y {mustBeNumeric, mustBeReal}
    k (1,1) double {mustBePositive}
    progressFcn = []
end

[nCells, nGenes] = size(X);
if numel(y) ~= nCells
    error('pkg:e_lrdetest:labelCount', ...
        ['Y has %d entries for %d cells. X is cells-by-genes, so Y needs ', ...
        'one label per row.'], numel(y), nCells);
end
y = double(y(:));
if numel(unique(y)) ~= 2
    error('pkg:e_lrdetest:needTwoGroups', ...
        'Y must take exactly two distinct values, found %d.', ...
        numel(unique(y)));
end

if issparse(X), X = full(X); end

% A logistic fit on a marker gene is often perfectly separable, which is
% the answer we want rather than a problem: the fit still maximises the
% likelihood, it just walks to the iteration limit complaining. Scoped, and
% restored on every path out.
warnState = warning();
restoreWarn = onCleanup(@() warning(warnState));
warning('off', 'stats:glmfit:IterationLimit');
warning('off', 'stats:glmfit:PerfectSeparation');
warning('off', 'stats:glmfit:IllConditioned');
warning('off', 'stats:glmfit:BadScaling');
warning('off', 'MATLAB:rankDeficientMatrix');

% Inf, not 0, for a gene that cannot be tested: MINK below takes the
% smallest, and an untested gene must never be picked as a marker.
pval = inf(nGenes, 1);

% A gene that does not vary across the cells carries no information about
% the grouping and makes the fit degenerate. Skipping them is also most of
% the running time on real data, where the majority of genes are undetected
% in any given selection.
testable = find(std(X, 0, 1) > 0);

reportEvery = max(1, floor(numel(testable)/100));
for t = 1:numel(testable)
    g = testable(t);
    modelData = table(X(:, g), y, 'VariableNames', {'GENE', 'Group'});
    try
        m1 = fitglm(modelData, 'Group ~ GENE', 'Distribution', 'binomial');
        m2 = fitglm(modelData, 'Group ~ 1', 'Distribution', 'binomial');
        stat = 2*(m1.LogLikelihood - m2.LogLikelihood);
        df = m1.NumPredictors - m2.NumPredictors;
        pval(g) = chi2cdf(stat, df, 'upper');
    catch
        % A fit that will not converge at all leaves this gene at Inf,
        % which excludes it rather than ranking it.
    end
    if ~isempty(progressFcn) && mod(t, reportEvery) == 0
        progressFcn(t/numel(testable));
    end
end

% Never ask for more than there is, and never return an untested gene.
k = min(floor(k), nGenes);
[~, idx] = mink(pval, k);
idx = idx(isfinite(pval(idx)));
idx = idx(:);

end
