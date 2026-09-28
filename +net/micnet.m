function A = micnet(X, opts)
%MICNET Gene network from pairwise maximal information coefficient.
%   A = net.micnet(X) returns the genes-by-genes MIC matrix for X, a
%   genes-by-cells expression matrix. MIC (Reshef et al., Science 2011)
%   scores any functional dependence, linear or not, on a 0-1 scale
%   meant to be comparable across relationship shapes.
%
%   The diagonal is left at 0, as in net.distcorrnet: A is an adjacency
%   matrix without self-loops.
%
%   Name-value arguments:
%     MaxCells - use at most this many cells, evenly spaced through the
%                columns of X (default 2000). One pair takes about 0.05 s
%                at 1000 cells and 0.2 s at 2000, and there are
%                G(G-1)/2 pairs for G genes.
%
%   Pairs are split across an already-open parallel pool, if there is
%   one; none is started here. Unrelated genes score about 0.1, not 0,
%   so threshold A relative to that floor rather than to zero.
%
%   See also run.ml_mine, net.distcorrnet, net.minet.

arguments
    X {mustBeNumeric}
    opts.MaxCells (1, 1) double {mustBeInteger, mustBePositive} = 2000
end

X = double(full(X));
[numGenes, numCells] = size(X);
keep = unique(round(linspace(1, numCells, min(numCells, opts.MaxCells))));
X = X(:, keep);

[rows, cols] = find(triu(true(numGenes), 1));
numPairs = numel(rows);
scores = zeros(numPairs, 1);
numWorkers = poolSize();
parfor (k = 1:numPairs, numWorkers)
    s = run.ml_mine(X(rows(k), :), X(cols(k), :));
    scores(k) = s.mic;
end

A = zeros(numGenes);
A(sub2ind([numGenes, numGenes], rows, cols)) = scores;
A = A + A.';
end

function n = poolSize()
% Workers in the current pool, or 0 -- which makes parfor run serially and
% never auto-create a pool. GCP is absent without Parallel Computing
% Toolbox, so its failure means no pool.
n = 0;
try
    p = gcp("nocreate");
    if ~isempty(p)
        n = p.NumWorkers;
    end
catch
    % No Parallel Computing Toolbox: run serially.
end
end
