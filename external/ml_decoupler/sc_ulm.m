function [scores, pvals, terms] = sc_ulm(sce, net, options)
% SC_ULM  Cell-level enrichment scoring via Univariate Linear Model (ULM).
%
% Equivalent to decoupleR / decoupler-py dc.mt.ulm(). For each cell and
% each gene set (term), fits:
%
%   expression_g = intercept + beta * weight_g + epsilon
%
% and returns the OLS t-statistic of beta as the enrichment score.
% Larger magnitude = more significant enrichment; sign indicates over-
% (positive) or under-representation (negative) of the gene set.
%
% USAGE:
%   [scores, pvals, terms] = sc_ulm(sce, net)
%   [scores, pvals, terms] = sc_ulm(sce, net, MinGenes=5, Verbose=true)
%
% INPUTS:
%   sce  - SingleCellExperiment object. sce.X must contain log-normalized
%          counts (genes x cells). Use sc_norm then log1p beforehand if needed.
%   net  - Table with columns:
%            source  (string) — gene set / TF / pathway name
%            target  (string) — gene symbol
%            weight  (double) — regulatory weight (optional; default = 1)
%          Typical sources: sc_net_progeny(), sc_net_panglaodb(),
%          sc_net_collectri(), or any custom gene-set table.
%
% OUTPUTS:
%   scores - n_cells x n_terms matrix of enrichment t-statistics
%   pvals  - n_cells x n_terms matrix of two-sided p-values (NaN for skipped terms)
%   terms  - 1 x n_terms string array of term names
%
% OPTIONS:
%   MinGenes (default 3) - terms with fewer overlapping genes are skipped
%   Verbose  (default true) - print progress summary
%
% EXAMPLE:
%   net = sc_net_progeny("human");
%   [scores, pvals, terms] = sc_ulm(sce, net, MinGenes=5);
%   T = sc_rankby_group(scores, terms, sce.c_cluster_id);
%
% NOTE ON WEIGHTS:
%   ULM measures co-variation of gene expression with gene weights within each
%   gene set — not overall expression level. Enrichment is detected when
%   highly-weighted genes are expressed proportionally more than low-weighted
%   genes. Uniform weights (all identical) make the model degenerate; use
%   weighted databases (PROGENy, PanglaoDB with sensitivity scores, DoRothEA).
%
% REFERENCE:
%   Badia-i-Mompel et al., Bioinformatics 2022.
%   https://doi.org/10.1093/bioinformatics/btac832

arguments
    sce SingleCellExperiment
    net table
    options.MinGenes (1,1) double {mustBePositive, mustBeInteger} = 3
    options.Verbose  (1,1) logical = true
end

% --- Validate network table ---
for col = ["source", "target"]
    if ~ismember(col, net.Properties.VariableNames)
        error('sc_ulm:missingColumn', ...
            'Network table must have a "%s" column.', col);
    end
end
if ~ismember("weight", net.Properties.VariableNames)
    net.weight = ones(height(net), 1);
end
net.source = string(net.source);
net.target = string(net.target);
net.weight = double(net.weight);

% --- Expression matrix (genes x cells), dense double ---
X = full(double(sce.X));
geneNames = upper(string(sce.g));
[n_genes, n_cells] = size(X);

if options.Verbose
    fprintf('SC_ULM: %d genes x %d cells\n', n_genes, n_cells);
end

% --- Log-normalize if data looks like raw counts (max per cell >> 20) ---
if max(X(:)) > 50
    warning('sc_ulm:unnormalized', ...
        ['sce.X contains large values (max=%.0f). ' ...
         'Apply sc_norm and log1p before calling sc_ulm.'], max(X(:)));
end

% --- Build term list ---
terms = unique(net.source, 'stable');
n_terms = numel(terms);

scores = nan(n_cells, n_terms);
pvals  = ones(n_cells, n_terms);
n_used = 0;

for k = 1:n_terms
    mask = net.source == terms(k);
    termTargets  = upper(net.target(mask));
    termWeights  = net.weight(mask);

    % Intersect with measured genes (case-insensitive)
    [~, iGene, iNet] = intersect(geneNames, termTargets, 'stable');
    n_overlap = numel(iGene);

    if n_overlap < options.MinGenes
        continue
    end

    X_t = X(iGene, :);           % n_overlap x n_cells
    w_t = termWeights(iNet);      % n_overlap x 1

    [t_stat, p_val] = i_ulm_ols(X_t, w_t);
    scores(:, k) = t_stat';
    pvals(:, k)  = p_val';
    n_used = n_used + 1;
end

if options.Verbose
    fprintf('SC_ULM: scored %d/%d terms (skipped %d with <%d genes).\n', ...
        n_used, n_terms, n_terms - n_used, options.MinGenes);
end

end


% -------------------------------------------------------------------------
function [t_scores, p_vals] = i_ulm_ols(X_t, w_t)
% Vectorized OLS ULM for one gene set across all cells simultaneously.
%
% Fits: X_t(:, i) = intercept + beta * w_t + epsilon  for each cell i.
% Returns t-statistic and p-value of beta (the enrichment coefficient).
%
% This removes the intercept by mean-centering, then applies closed-form OLS:
%   beta  = (w_c' * X_c) / (w_c' * w_c)
%   t     = beta / (sqrt(SSE/(n-2)) / sqrt(w_c' * w_c))

n_genes = size(X_t, 1);

% Mean-center each cell across its genes; center the weight vector
X_c = X_t - mean(X_t, 1);   % n_genes x n_cells
w_c = w_t - mean(w_t);       % n_genes x 1

ww = sum(w_c .^ 2);

if ww < eps
    t_scores = zeros(1, size(X_t, 2));
    p_vals   = ones(1, size(X_t, 2));
    return
end

% OLS coefficients for all cells at once
beta = (w_c' * X_c) / ww;          % 1 x n_cells

% Residuals and sum-of-squared errors per cell
residuals = X_c - w_c * beta;      % n_genes x n_cells
SSE = sum(residuals .^ 2, 1);      % 1 x n_cells

% Standard error of beta; guard against zero SE
se = sqrt(SSE / max(n_genes - 2, 1)) / sqrt(ww);
se = max(se, eps);

% t-statistics and two-sided p-values
t_scores = beta ./ se;
p_vals   = 2 * (1 - tcdf(abs(t_scores), max(n_genes - 2, 1)));

end
