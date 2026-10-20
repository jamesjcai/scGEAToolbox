function model = e_genecorrmodel(X, g, opts)
%E_GENECORRMODEL  What PKG.E_GENECORRROW needs to correlate one gene with every gene.
%
%   model = pkg.e_genecorrmodel(X, g) takes the raw count matrix X
%   (genes x cells) and its gene names G and returns a struct with fields
%
%     Genes     the gene names, a string column
%     Mean      mean log1p library-size-normalised expression of each gene
%     Fraction  fraction of cells in which each gene is detected
%     Pearson   struct with the matrix the Pearson correlation is taken on
%               (log1p normalised expression), its row means and row
%               standard deviations
%     Spearman  the same for the Spearman correlation, whose matrix holds
%               each gene's ranks across cells; [] unless asked for
%
%   model = pkg.e_genecorrmodel(X, g, Spearman=true) also ranks the genes.
%   That doubles the memory, so it is off by default.
%
%   Both matrices stay sparse. Ranking a sparse gene does not need a dense
%   copy: its zeros tie, all at rank (z+1)/2 for z zeros, and its nonzero
%   values rank above them, from z+1 up. Subtracting (z+1)/2 from every rank
%   leaves the zeros at zero and does not change any correlation, so the
%   ranks are stored that way.
%
%   See also PKG.E_GENECORRROW, GUI.CALLBACK_GENECORRPARTNERS.

arguments
    X {mustBeNumeric}
    g
    opts.Spearman (1, 1) logical = false
end

[numGenes, numCells] = size(X);
if numel(g) ~= numGenes
    error("e_genecorrmodel:SizeMismatch", ...
        "G has %d names but X has %d genes. Pass one name per row of X.", ...
        numel(g), numGenes);
end

Y = sparse(double(log1p(sc_norm(X))));
Y(isnan(Y)) = 0;

counts = full(sum(Y ~= 0, 2));
spearman = [];
if opts.Spearman
    spearman = i_moments(i_ranks(Y, counts));
end

model = struct("Genes", string(g(:)), ...
    "Mean", full(mean(Y, 2)), ...
    "Fraction", counts/numCells, ...
    "Pearson", i_moments(Y), ...
    "Spearman", spearman);
end

function R = i_ranks(Y, counts)
% Rank each gene's nonzero values among themselves, ties averaged, in the
% shifted form described above: rank + (z - 1)/2 for a gene with z zeros.
[numGenes, numCells] = size(Y);
[cellIdx, geneIdx, v] = find(Y.');
ranks = zeros(size(v));
stop = cumsum(counts);
start = stop - counts + 1;
for k = find(counts > 0).'
    span = start(k):stop(k);
    numZeros = numCells - counts(k);
    ranks(span) = tiedrank(v(span)) + (numZeros - 1)/2;
end
R = sparse(geneIdx, cellIdx, ranks, numGenes, numCells);
end

function m = i_moments(M)
% The matrix with its row means and sample standard deviations.
n = size(M, 2);
mu = full(sum(M, 2))/n;
sumSquares = full(sum(M.^2, 2));
sd = sqrt(max(sumSquares - n*mu.^2, 0)/(n - 1));
m = struct("Matrix", M, "Mean", mu, "Std", sd);
end
