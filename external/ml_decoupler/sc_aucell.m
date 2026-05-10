function [scores, terms] = sc_aucell(sce, net, options)
% SC_AUCELL  Cell-level gene-set enrichment via Area Under the recovery Curve.
%
% For each cell, genes are ranked by expression (rank 1 = highest). For each
% gene set, the AUC measures how many set members fall within the top-ranked
% genes versus a random expectation. AUC = 1 means all set genes are at the
% very top; AUC ~ 0 means no enrichment.
%
% AUCell does NOT require weighted gene sets — only membership matters.
% Equivalent to decoupleR / decoupler-py dc.mt.aucell().
%
% USAGE:
%   [scores, terms] = sc_aucell(sce, net)
%   [scores, terms] = sc_aucell(sce, net, TopGenes=500, MinGenes=5)
%
% INPUTS:
%   sce  - SingleCellExperiment with log-normalized or raw counts in sce.X
%   net  - Table with columns: source (term), target (gene), weight
%          Weights are ignored — only membership is used.
%
% OUTPUTS:
%   scores - n_cells x n_terms matrix of AUC values in [0, 1]
%   terms  - 1 x n_terms string array of term names
%
% OPTIONS:
%   TopGenes (default: ceil(0.05 * n_genes)) - number of top genes per cell
%            to consider as the "expressed" pool (5% rule from Aibar et al.)
%   MinGenes (default 3) - terms with fewer overlapping genes are skipped
%
% ALGORITHM:
%   For each cell and gene set with k members in the top TopGenes genes:
%     AUC  = TopGenes * k - sum(positions)
%     max  = TopGenes * k - k*(k+1)/2         (best case: genes at rank 1..k)
%     score = AUC / max
%   This is a closed-form equivalent of the running-integral computation.
%
% REFERENCE:
%   Aibar et al., Nature Methods 2017.
%   https://doi.org/10.1038/nmeth.4463

arguments
    sce SingleCellExperiment
    net table
    options.TopGenes (1,1) double = 0      % 0 = auto (5% of genes)
    options.MinGenes (1,1) double {mustBePositive, mustBeInteger} = 3
end

X = full(double(sce.X));
[n_genes, n_cells] = size(X);
geneNames = upper(string(sce.g));

n_up = options.TopGenes;
if n_up <= 0
    n_up = max(2, ceil(0.05 * n_genes));
end
n_up = min(n_up, n_genes);

% Validate and normalise network
for col = ["source", "target"]
    assert(ismember(col, net.Properties.VariableNames), ...
        'sc_aucell: network table must have a "%s" column.', col);
end
net.source = string(net.source);
net.target = upper(string(net.target));

terms_all = unique(net.source, 'stable');
n_terms   = numel(terms_all);

% --- Precompute global rank matrix (rank 1 = highest expression per cell) ---
% ranks(g, c) = rank of gene g in cell c (integer, 1-based)
[~, sortIdx] = sort(X, 1, 'descend');
ranks = zeros(n_genes, n_cells, 'uint32');
colSub = repmat(uint32(1:n_cells), n_genes, 1);
ranks(sub2ind([n_genes, n_cells], sortIdx, colSub)) = ...
    repmat(uint32(1:n_genes)', 1, n_cells);

scores    = zeros(n_cells, n_terms);
n_kept    = 0;
keepTerm  = false(n_terms, 1);

for k = 1:n_terms
    mask = net.source == terms_all(k);
    tgts = net.target(mask);

    [~, iGene, ~] = intersect(geneNames, tgts, 'stable');
    if numel(iGene) < options.MinGenes
        continue
    end
    keepTerm(k) = true;
    n_kept = n_kept + 1;

    R_k = double(ranks(iGene, :));   % n_overlap x n_cells

    % Number of set genes in top n_up per cell
    in_top  = R_k <= n_up;           % logical n_overlap x n_cells
    k_c     = sum(in_top, 1);        % 1 x n_cells

    % AUC = n_up * k_c - sum of ranks for in-top genes
    sum_pos = sum(R_k .* in_top, 1); % 1 x n_cells

    auc     = double(n_up) * k_c - sum_pos;
    max_auc = double(n_up) * k_c - k_c .* (k_c + 1) / 2;

    score = zeros(1, n_cells);
    valid = max_auc > 0;
    score(valid) = auc(valid) ./ max_auc(valid);
    scores(:, k) = score';
end

scores = scores(:, keepTerm);
terms  = terms_all(keepTerm)';

fprintf('SC_AUCELL: scored %d/%d terms (TopGenes=%d, skipped %d with <%d genes).\n', ...
    n_kept, n_terms, n_up, n_terms - n_kept, options.MinGenes);

end
