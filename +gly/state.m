function [cs, setnames, ncommon] = state(X, genelist, methodid, minGenes)
%GLY.STATE  Per-cell glycobiological state from a glycogene set collection.
%
%   [cs, setnames] = GLY.STATE(X, genelist) scores every cell against
%   each gene set in the curated glycobiology collection (GLY.GENESETS)
%   and returns an M-by-N matrix CS of module scores (M glyco modules x
%   N cells). Together, a cell's row of scores defines its glycobiological
%   state - the relative activity of N-/O-glycosylation, sialylation,
%   fucosylation, glycosaminoglycan and glycosphingolipid metabolism, glycan
%   degradation, and glycan-recognition (lectin) programs.
%
%   USAGE:
%     [cs, setnames] = gly.state(sce.X, sce.g);
%     [cs, setnames] = gly.state(X, genelist, methodid);
%
%   INPUTS:
%     X        - genes-by-cells expression matrix (raw or normalized counts)
%     genelist - G-by-1 gene symbols (length = rows of X)
%     methodid - scoring method forwarded to SC_CELLSCORE (default 2,
%                AddModuleScore): 1 = UCell, 2 = AddModuleScore, 3 = AUCell
%     minGenes - minimum number of a set's genes that must be present in the
%                data for the set to be scored (default 3). Sets below the
%                threshold are dropped from the output.
%
%   OUTPUTS:
%     cs       - M-by-N glyco-module score matrix (rows follow setnames)
%     setnames - M-by-1 string array of the scored module names
%     ncommon  - M-by-1 count of each module's genes found in the data
%
% see also: GLY.GENESETS, SC_CELLSCORE, SC_PATHWAYACTIVITY

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    methodid (1, 1) double = 2
    minGenes (1, 1) double = 3
end

if numel(genelist) ~= size(X, 1)
    error("GENELIST length (%d) must equal the number of rows of X (%d).", ...
        numel(genelist), size(X, 1));
end

[setmatrx, allnames, setgenes] = gly.genesets();

nSets = numel(allnames);
nCells = size(X, 2);
csAll = NaN(nSets, nCells);
ncommonAll = zeros(nSets, 1);

for k = 1:nSets
    tgsPos = setgenes(setmatrx(k, :));
    ncommonAll(k) = sum(ismember(upper(genelist), upper(tgsPos)));
    if ncommonAll(k) < minGenes
        continue;
    end
    csAll(k, :) = sc_cellscore(X, genelist, tgsPos, [], methodid);
end

keep = ncommonAll >= minGenes;
cs = csAll(keep, :);
setnames = allnames(keep);
ncommon = ncommonAll(keep);

if isempty(setnames)
    warning("No glyco module had at least %d genes present in the data.", minGenes);
end

end
