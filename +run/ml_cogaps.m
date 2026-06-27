function [A, P, Tmarkers] = ml_cogaps(X, g, nPatterns, nIterations, sparseOpt)
% ML_COGAPS - MATLAB-native NMF pattern discovery (CoGAPS-like)
%
% A lightweight MATLAB reimplementation of the CoGAPS *workflow* using the
% built-in NNMF function. It is *not* the Bayesian MCMC algorithm of CoGAPS:
% it computes a classical (Frobenius-norm) non-negative factorization and the
% PatternMarker gene ranking on top of it. Use this when R/Bioconductor is
% not available; use RUN.R_COGAPS for the true Bayesian method.
%
% [A, P, Tmarkers] = run.ml_cogaps(X, g, nPatterns, nIterations, sparseOpt)
%   X           - gene-by-cell expression matrix (counts; nonnegative)
%   g           - gene names (string array, length = size(X,1))
%   nPatterns   - number of patterns to learn (default 8)
%   nIterations - max NNMF iterations (default 1000)
%   sparseOpt   - prefer sparser factors (multiplicative update) (default true)
% Returns:
%   A        - genes-by-patterns feature loadings (A matrix)
%   P        - cells-by-patterns sample factors (P matrix)
%   Tmarkers - PatternMarker table (gene, pattern, rank, score)
%
% See also: SC_NMFPATTERN, RUN.R_COGAPS, NNMF

arguments
    X {mustBeNumeric}
    g = []
    nPatterns (1, 1) double = 8
    nIterations (1, 1) double = 1000
    sparseOpt (1, 1) logical = true
end

if issparse(X), X = full(X); end
X = double(X);
if isempty(g), g = pkg.i_defaultgenenames(size(X, 1), "G"); end
g = string(g(:));
if numel(g) ~= size(X, 1)
    error('ml_cogaps:DimensionMismatch', ...
        'Length of g (%d) must match number of rows in X (%d).', ...
        numel(g), size(X, 1));
end

% NNMF requires nonnegative input. Normalize and log-transform counts for
% well-behaved patterns (standard practice for single-cell NMF).
X = max(X, 0);
Xin = log1p(sc_norm(X));

% nPatterns must be smaller than both dimensions of the matrix.
nPatterns = min(nPatterns, min(size(Xin)) - 1);
if nPatterns < 2
    error('ml_cogaps:TooFewPatterns', ...
        'Not enough genes/cells to learn at least 2 patterns.');
end

if sparseOpt
    algo = 'mult';   % multiplicative update tends to give sparser factors
else
    algo = 'als';    % alternating least squares is faster
end
opt = statset('MaxIter', nIterations, 'Display', 'off');
[W, H] = nnmf(Xin, nPatterns, 'replicates', 5, 'algorithm', algo, 'options', opt);

A = W;        % genes x patterns
P = H';       % cells x patterns

Tmarkers = i_patternmarkers(A, g);
end

function T = i_patternmarkers(A, g)
% PatternMarker statistic: assign each gene to the pattern it is most
% specific to, then rank genes within each pattern (CoGAPS threshold="all").
k = size(A, 2);

% Scale each pattern to [0,1], then normalize each gene across patterns.
An = A ./ max(max(A, [], 1), eps);
rs = sum(An, 2);
rs(rs == 0) = eps;
An = An ./ rs;

% Distance of each gene to the "perfect marker" reference of each pattern.
dist = zeros(size(A, 1), k);
for p = 1:k
    ref = zeros(1, k);
    ref(p) = 1;
    dist(:, p) = sqrt(sum((An - ref).^2, 2));
end

[score, bestPattern] = min(dist, [], 2);
gene = g(:);
pattern = bestPattern;
T = table(gene, pattern, score);
T = sortrows(T, {'pattern', 'score'});

rank = zeros(height(T), 1);
for p = 1:k
    ix = T.pattern == p;
    rank(ix) = 1:nnz(ix);
end
T.rank = rank;
T = T(:, {'gene', 'pattern', 'rank', 'score'});
end
