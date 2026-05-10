function [scores, pvals, terms] = sc_zscore(sce, net, options)
% SC_ZSCORE  Cell-level enrichment scoring via weighted z-score normalization.
%
% Computes the enrichment of a gene set as the (weighted) mean expression
% of its target genes, scaled by the within-cell standard deviation and
% normalized by the square root of the number of target genes.
%
% Two flavors are available:
%   RoKAI: ES = ((µ_set - µ_global) × √m) / σ_global
%   KSEA:  ES = (µ_set × √m) / σ_global
%
% Equivalent to decoupler-py dc.mt.zscore().
%
% USAGE:
%   [scores, pvals, terms] = sc_zscore(sce, net)
%   [scores, pvals, terms] = sc_zscore(sce, net, Flavor="KSEA", MinGenes=3)
%
% INPUTS:
%   sce  - SingleCellExperiment with log-normalized counts in sce.X
%   net  - Table with columns: source (term), target (gene), weight
%
% OUTPUTS:
%   scores - n_cells x n_terms matrix of z-scores
%   pvals  - n_cells x n_terms matrix of two-sided p-values (standard normal)
%   terms  - 1 x n_terms string array of term names
%
% OPTIONS:
%   Flavor   ("RoKAI" | "KSEA", default "RoKAI")
%            RoKAI: subtracts per-cell global mean before scaling
%            KSEA:  uses raw weighted mean without global-mean correction
%   MinGenes (default 3) - skip terms with fewer overlapping genes
%
% ALGORITHM:
%   1. Build weighted adjacency matrix A (n_genes x n_terms)
%   2. Per cell c:
%      σ_c   = std(X[:, c])
%      µ_c   = mean(X[:, c])      [RoKAI only]
%      µ_set = (A' * X[:, c]) / sum(|A|)    [weighted mean per term]
%      m     = sqrt(|{genes in set}|)
%      ES    = (µ_set - µ_c) * m / σ_c      [RoKAI]
%           or µ_set * m / σ_c              [KSEA]
%   3. p-value from standard normal: 2 * (1 - Φ(|ES|))
%
% REFERENCE:
%   Hernandez-Armenta et al., Bioinformatics 2017 (KSEA).
%   https://doi.org/10.1093/bioinformatics/btx082
%   Yılmaz et al., Cell Systems 2021 (RoKAI).
%   https://doi.org/10.1016/j.cels.2021.08.012

arguments
    sce SingleCellExperiment
    net table
    options.Flavor   (1,1) string  = "RoKAI"
    options.MinGenes (1,1) double  {mustBePositive, mustBeInteger} = 3
end

flavor = lower(options.Flavor);
assert(ismember(flavor, ["rokai", "ksea"]), ...
    'sc_zscore: Flavor must be "RoKAI" or "KSEA". Got: %s', options.Flavor);

X = full(double(sce.X));
[~, n_cells] = size(X);
geneNames = string(sce.g);

for col = ["source", "target"]
    assert(ismember(col, net.Properties.VariableNames), ...
        'sc_zscore: network table must have a "%s" column.', col);
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
    fprintf('SC_ZSCORE: no terms passed MinGenes=%d filter.\n', options.MinGenes);
    return
end

% --- Per-cell statistics ---
mu_global  = mean(X, 1);    % 1 x n_cells
sig_global = std(X, 0, 1);  % 1 x n_cells (ddof=1 equivalent)
sig_global = max(sig_global, eps);

% --- Weighted mean of targets per term per cell ---
% adj is sparse n_genes x n_terms
% X' * adj  =>  n_cells x n_terms  (sum of expression * weight per term per cell)
% divide by sum(|adj|, 1)  =>  normalize to weighted mean
w_sum  = full(sum(abs(adj), 1));        % 1 x n_terms
w_sum  = max(w_sum, eps);
mu_set = (X' * full(adj)) ./ w_sum;    % n_cells x n_terms

% --- Number of non-zero genes per term ---
n_k = sqrt(full(sum(adj ~= 0, 1)));    % 1 x n_terms (sqrt of set size)

% --- Enrichment score ---
if strcmp(flavor, "rokai")
    % Subtract per-cell global mean
    scores = ((mu_set - mu_global') .* n_k) ./ sig_global';   % n_cells x n_terms
else
    % KSEA: no mean subtraction
    scores = (mu_set .* n_k) ./ sig_global';
end

% --- Two-sided p-values from standard normal ---
pvals = 2 * (1 - normcdf(abs(scores)));

fprintf('SC_ZSCORE (%s): scored %d terms across %d cells.\n', ...
    options.Flavor, n_terms, n_cells);

end
