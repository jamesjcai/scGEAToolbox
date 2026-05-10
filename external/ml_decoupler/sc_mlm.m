function [scores, pvals, terms] = sc_mlm(sce, net, options)
% SC_MLM  Cell-level enrichment via Multivariate Linear Model (MLM).
%
% Fits a single multivariate OLS regression per cell where all gene-set
% weights enter simultaneously as predictors:
%
%   expression_g = β_0 + β_1·w_{g,1} + β_2·w_{g,2} + … + β_T·w_{g,T} + ε
%
% The t-statistic of each β_k is the enrichment score for gene set k.
% Unlike ULM (which fits each gene set in isolation), MLM accounts for
% shared genes across sets, making it better suited to dense networks like
% DoRothEA where many TFs share target genes.
%
% Equivalent to decoupler-py dc.mt.mlm().
%
% USAGE:
%   [scores, pvals, terms] = sc_mlm(sce, net)
%   [scores, pvals, terms] = sc_mlm(sce, net, MinGenes=5)
%
% INPUTS:
%   sce  - SingleCellExperiment with log-normalized counts in sce.X
%   net  - Table with columns: source (term), target (gene), weight
%
% OUTPUTS:
%   scores - n_cells x n_terms matrix of t-statistics (enrichment scores)
%   pvals  - n_cells x n_terms matrix of two-sided p-values
%   terms  - 1 x n_terms string array of term names
%
% OPTIONS:
%   MinGenes (default 3) - skip terms with fewer overlapping genes
%   ReturnCoef (false)   - if true, return OLS coefficients instead of t-stats
%
% ALGORITHM:
%   1. Build design matrix A = [1, adj]  (n_genes × 1+n_terms)
%      where adj(g, k) = weight of gene g in set k (0 if not a member)
%   2. Precompute (A'A)⁻¹ (same for all cells — computed once)
%   3. For all cells simultaneously: coef = (A'A)⁻¹ · A' · X
%      (n_terms+1 × n_cells via a single matrix multiply)
%   4. Per cell: residuals = X - A·coef; MSE = SSE / df
%   5. t-stat for term k = coef_k / (√MSE · √[(A'A)⁻¹]_{kk})
%
% NOTES:
%   - Requires n_genes >> n_terms for well-defined degrees of freedom.
%     If df ≤ 0, an error is raised.
%   - When many gene sets share few members, MLM may produce inflated
%     scores relative to ULM; benchmark both for your dataset.
%
% REFERENCE:
%   Badia-i-Mompel et al., Bioinformatics 2022.
%   https://doi.org/10.1093/bioinformatics/btac832

arguments
    sce SingleCellExperiment
    net table
    options.MinGenes   (1,1) double {mustBePositive, mustBeInteger} = 3
    options.ReturnCoef (1,1) logical = false
end

X = full(double(sce.X));
[n_genes, n_cells] = size(X);
geneNames = string(sce.g);

for col = ["source", "target"]
    assert(ismember(col, net.Properties.VariableNames), ...
        'sc_mlm: network table must have a "%s" column.', col);
end
if ~ismember("weight", net.Properties.VariableNames)
    net.weight = ones(height(net), 1);
end

% Build adjacency matrix (n_genes x n_terms)
[adj, terms] = i_build_adj(geneNames, net, options.MinGenes);
n_terms = numel(terms);

if n_terms == 0
    scores = zeros(n_cells, 0);
    pvals  = ones(n_cells, 0);
    fprintf('SC_MLM: no terms passed MinGenes=%d filter.\n', options.MinGenes);
    return
end

% Degrees of freedom
df = n_genes - n_terms - 1;
if df <= 0
    error('sc_mlm:insufficientDOF', ...
        ['Degrees of freedom = %d (n_genes=%d, n_terms=%d). ' ...
         'Need n_genes >> n_terms. Reduce MinGenes or use sc_ulm instead.'], ...
        df, n_genes, n_terms);
end

adj_full = full(adj);

% Design matrix: prepend intercept column
A = [ones(n_genes, 1), adj_full];   % n_genes x (1+n_terms)

% Precompute (A'A)^{-1} once — shape (1+n_terms) x (1+n_terms)
AtA     = A' * A;
invAtA  = inv(AtA);

% OLS coefficients for ALL cells: (1+n_terms) x n_cells
coef = AtA \ (A' * X);

% Residuals and per-cell MSE
residuals = X - A * coef;           % n_genes x n_cells
SSE       = sum(residuals .^ 2, 1); % 1 x n_cells
mse       = SSE / df;               % 1 x n_cells

if options.ReturnCoef
    scores = coef(2:end, :)';       % n_cells x n_terms
    pvals  = ones(n_cells, n_terms);
    fprintf('SC_MLM: returning OLS coefficients (not t-stats) for %d terms.\n', n_terms);
    return
end

% Diagonal of (A'A)^{-1} for term blocks only (skip intercept row/col 1)
diag_inv = diag(invAtA(2:end, 2:end));   % n_terms x 1

% Standard errors: se(k, c) = sqrt(mse(c) * diag_inv(k))
% mse is 1 x n_cells, diag_inv is n_terms x 1 → broadcasting → n_terms x n_cells
se = sqrt(diag_inv .* mse);          % n_terms x n_cells
se = max(se, eps);

% t-statistics and two-sided p-values
t_stats = coef(2:end, :) ./ se;     % n_terms x n_cells
scores  = t_stats';                  % n_cells x n_terms
pvals   = 2 * (1 - tcdf(abs(scores), df));

fprintf('SC_MLM: scored %d terms across %d cells (df=%d).\n', ...
    n_terms, n_cells, df);

end
