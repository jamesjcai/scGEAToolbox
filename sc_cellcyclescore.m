function [ScoreV, T] = sc_cellcyclescore(X, g, options)
% Score cell cycle phases
%  [phase, T] = SC_CELLCYCLESCORE(X, g) scores each cell for the S and G2M
%  gene sets (Seurat AddModuleScore) and calls G1 when both scores are <= 0,
%  otherwise the phase with the higher score.
%
%  [...] = SC_CELLCYCLESCORE(X, g, ReferenceX=Xref, ReferenceGenes=gref)
%  takes the scoring baseline from a reference sample, such as untreated
%  control cells, instead of from X. The baseline is relative to the cells
%  it is built from, so a sample where most cells sit in one phase -- a
%  drug-arrested population, or a subset picked by an earlier call -- moves
%  its own baseline and is miscalled. See SC_CELLSCORE.

arguments
    X
    g
    options.ReferenceX = []
    options.ReferenceGenes = []
end

% Define path to cell cycle gene list
pw1 = fileparts(mfilename('fullpath'));
wrkpth = fullfile(pw1, 'assets', 'CellScores', 'cellcyclegenes.xlsx');

% Read gene table from the file
T = readtable(wrkpth, ...
    'ReadVariableNames', true, 'FileType', 'spreadsheet', ...
    'Sheet', 'Regev_cell_cycle_genes');

% Extract S and G2M phase genes
sgenes = string(T.S);
sgenes = sgenes(strlength(sgenes) > 0);
g2mgenes = string(T.G2M);
g2mgenes = g2mgenes(strlength(g2mgenes) > 0);

% Calculate scores for S and G2M phases
refArgs = {"ReferenceX", options.ReferenceX, "ReferenceGenes", options.ReferenceGenes};
score_S = sc_cellscore(X, g, sgenes, [], 2, refArgs{:});
score_G2M = sc_cellscore(X, g, g2mgenes, [], 2, refArgs{:});

% Assign cell cycle phase based on scores
if all(isnan(score_S)) || all(isnan(score_G2M))
    ScoreV = string(repmat('unknown', size(X, 2), 1));
else
    ScoreV = string(repmat('G1', size(X, 2), 1));
    C = [score_S, score_G2M];
    [~, I] = max(C, [], 2);
    Cx = C > 0;
    i = sum(Cx, 2) > 0;
    ScoreV(i & I == 1) = "S";
    ScoreV(i & I == 2) = "G2M";
end

% Return a table if requested
if nargout > 1
    T = table(score_S, score_G2M, ScoreV);
end
end
