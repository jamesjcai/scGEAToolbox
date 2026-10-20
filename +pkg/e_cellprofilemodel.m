function model = e_cellprofilemodel(X, g, opts)
%E_CELLPROFILEMODEL  What PKG.E_CELLPROFILEROW needs to profile one cell's genes.
%
%   model = pkg.e_cellprofilemodel(X, g) takes the raw count matrix X
%   (genes x cells) and its gene names G and returns a struct with fields
%
%     Genes       the gene names, a string column
%     Expression  genes x cells sparse log1p library-size-normalised
%                 expression
%     Mean        each gene's mean expression over all cells
%     Fraction    fraction of cells in which each gene is detected
%     Searcher    a KDTREESEARCHER over the cells' principal-component
%                 scores (PKG.E_CELLPCS), to find a cell's neighbours
%
%   model = pkg.e_cellprofilemodel(..., NumComponents=NPC) sets the number
%   of components (default 30).
%
%   See also PKG.E_CELLPROFILEROW, GUI.CALLBACK_CELLGENEPROFILE.

arguments
    X {mustBeNumeric}
    g
    opts.NumComponents (1, 1) double {mustBePositive, mustBeInteger} = 30
end

[numGenes, numCells] = size(X);
if numel(g) ~= numGenes
    error("e_cellprofilemodel:SizeMismatch", ...
        "G has %d names but X has %d genes. Pass one name per row of X.", ...
        numel(g), numGenes);
end

Y = sparse(double(log1p(sc_norm(X))));
Y(isnan(Y)) = 0;
score = pkg.e_cellpcs(X, g, opts.NumComponents);

model = struct("Genes", string(g(:)), ...
    "Expression", Y, ...
    "Mean", full(mean(Y, 2)), ...
    "Fraction", full(sum(Y ~= 0, 2))/numCells, ...
    "Searcher", KDTreeSearcher(score));
end
