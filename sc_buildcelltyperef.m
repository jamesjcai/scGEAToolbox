function ref = sc_buildcelltyperef(X, g, labels, options)
%SC_BUILDCELLTYPEREF Summarise a labelled dataset into a cell-type reference.
%
%   ref = SC_BUILDCELLTYPEREF(X, g, labels) reduces an annotated count
%   matrix to the summary statistics SC_CELLTYPEANNOREF compares a query
%   against: one mean expression profile and one marker set per cell type.
%   No cell-level data is kept, so the result is small, fast to compare
%   against, and can be shared where the raw data cannot.
%
%   X       genes-by-cells raw counts
%   g       gene names, one per row of X
%   labels  cell type of each cell; "" and missing are left out
%
%   Name-value options:
%     Species    recorded in ref.Species ("")
%     Source     free text recorded in ref.Source, e.g. an atlas and version
%     MinCells   cell types with fewer cells are dropped (10)
%     NumMarkers, MinLog2FC, MinPctIn, MaxPctOut, Alpha
%                marker filter, passed to PKG.I_CLUSTERMARKERS
%
%   The reference is a struct:
%     Genes     G-by-1 string
%     Types     1-by-T string
%     Profile   G-by-T mean of log1p(library-size-normalised counts)
%     Markers   1-by-T cell of string columns
%     NumCells  1-by-T
%     Species, Source
%
%   A reference published only as tables -- a per-type average expression
%   matrix, a marker list, or both -- can be put in the same struct by hand.
%   Either Profile or Markers may be left empty; SC_CELLTYPEANNOREF then
%   uses whichever test the reference supports.
%
%   Example:
%     ref = sc_buildcelltyperef(atlas.X, atlas.g, atlas.c_cell_type_tx, ...
%         Species="mouse", Source="my atlas v2");
%     [lab, T] = sc_celltypeannoref(sce.X, sce.g, sce.c_cluster_id, ref);
%
%   See also SC_CELLTYPEANNOREF, SC_ANNOTATECELLS.

arguments
    X {mustBeNumeric}
    g {mustBeText}
    labels {mustBeText}
    options.Species (1, 1) string = ""
    options.Source (1, 1) string = ""
    options.MinCells (1, 1) double {mustBePositive} = 10
    options.NumMarkers (1, 1) double {mustBePositive} = 100
    options.MinLog2FC (1, 1) double = log2(1.5)
    options.MinPctIn (1, 1) double = 0.5
    options.MaxPctOut (1, 1) double = 0.25
    options.Alpha (1, 1) double = 0.05
end

g = string(g(:));
labels = string(labels(:));
if numel(g) ~= size(X, 1)
    error("sc_buildcelltyperef:GeneCount", ...
        "g has %d entries for %d rows of X. Pass one gene name per row.", ...
        numel(g), size(X, 1));
end
if numel(labels) ~= size(X, 2)
    error("sc_buildcelltyperef:LabelCount", ...
        "labels has %d entries for %d cells. Pass one label per cell.", ...
        numel(labels), size(X, 2));
end

isLabelled = ~ismissing(labels) & strlength(labels) > 0;
[c, types] = findgroups(labels(isLabelled));
counts = accumarray(c, 1)';
keep = counts >= options.MinCells;
if nnz(keep) < 2
    error("sc_buildcelltyperef:TooFewTypes", ...
        ['Fewer than two cell types have at least %d cells. Lower ', ...
         'MinCells or use a reference with more types.'], options.MinCells);
end
cells = find(isLabelled);
isKeptCell = keep(c);
cells = cells(isKeptCell);
[c, types] = findgroups(types(c(isKeptCell)));

X = X(:, cells);
Xn = log1p(sc_norm(X));
member = sparse(1:numel(c), c, 1, numel(c), numel(types));
numCells = full(sum(member, 1));
profile = full(Xn*member)./numCells;

markers = pkg.i_clustermarkers(X, g, c, NumMarkers=options.NumMarkers, ...
    MinLog2FC=options.MinLog2FC, MinPctIn=options.MinPctIn, ...
    MaxPctOut=options.MaxPctOut, Alpha=options.Alpha);

ref = struct("Genes", g, "Types", string(types(:))', "Profile", profile, ...
    "Markers", {markers}, "NumCells", numCells, ...
    "Species", options.Species, "Source", options.Source);
end
