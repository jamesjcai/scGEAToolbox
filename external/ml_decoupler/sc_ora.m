function [scores, pvals, terms] = sc_ora(sce, net, options)
% SC_ORA  Over-Representation Analysis per cell via Fisher's exact test.
%
% For each cell, the top TopGenes expressed genes are treated as the
% "significant" feature set. For each gene set (term), a 2x2 contingency
% table is constructed and a one-sided Fisher's exact test (hypergeometric)
% tests whether the set is over-represented in the top genes. The score is
% the log odds ratio (Haldane-Anscombe corrected).
%
% Equivalent to decoupler-py dc.mt.ora().
%
% USAGE:
%   [scores, pvals, terms] = sc_ora(sce, net)
%   [scores, pvals, terms] = sc_ora(sce, net, TopGenes=300, Background=15000)
%
% INPUTS:
%   sce  - SingleCellExperiment with log-normalized counts in sce.X
%   net  - Table with columns: source (term), target (gene), weight
%          Weights are ignored — only membership is used.
%
% OUTPUTS:
%   scores - n_cells x n_terms matrix of log odds ratios (Haldane corrected)
%            Positive = over-represented, negative = under-represented
%   pvals  - n_cells x n_terms matrix of p-values (right-sided Fisher's test)
%   terms  - 1 x n_terms string array of term names
%
% OPTIONS:
%   TopGenes   (default ceil(0.05 * n_genes)) - top expressed genes per cell
%   Background (default 0) - background gene pool size; 0 = use n_genes
%   MinGenes   (default 3) - skip terms with fewer overlapping genes
%
% CONTINGENCY TABLE (per cell c, gene set F):
%
%               In set F | Not in F
%   Top genes      a     |    b       (a + b = TopGenes)
%   Other genes    c     |    d
%   Total      |F ∩ bg|  | |bg| - |F|  (rows sum to Background)
%
%   log OR = log((a+0.5)(d+0.5) / ((b+0.5)(c+0.5)))
%   p      = 1 - hygecdf(a-1, Background, |F|, TopGenes)
%
% REFERENCE:
%   Subramanian et al., PNAS 2005 (ORA framework).
%   Fisher (1922); Haldane-Anscombe correction for log OR.

arguments
    sce SingleCellExperiment
    net table
    options.TopGenes   (1,1) double = 0     % 0 = auto (5% of genes)
    options.Background (1,1) double = 0     % 0 = use n_genes
    options.MinGenes   (1,1) double {mustBePositive, mustBeInteger} = 3
end

X = full(double(sce.X));
[n_genes, n_cells] = size(X);
geneNames = upper(string(sce.g));

n_up = options.TopGenes;
if n_up <= 0
    n_up = max(2, ceil(0.05 * n_genes));
end
n_up = min(n_up, n_genes);

n_bg = options.Background;
if n_bg <= 0
    n_bg = n_genes;
end
n_bg = max(n_bg, n_up);

for col = ["source", "target"]
    assert(ismember(col, net.Properties.VariableNames), ...
        'sc_ora: network table must have a "%s" column.', col);
end
net.source = string(net.source);
net.target = upper(string(net.target));

terms_all = unique(net.source, 'stable');
n_terms   = numel(terms_all);

% --- Precompute top-gene indicator matrix ---
% in_top(g, c) = true if gene g is in the top n_up genes of cell c
[~, sortIdx] = sort(X, 1, 'descend');
top_idx = sortIdx(1:n_up, :);          % n_up x n_cells (gene indices)
in_top  = false(n_genes, n_cells);
for c = 1:n_cells
    in_top(top_idx(:, c), c) = true;
end

scores   = zeros(n_cells, n_terms);
pvals    = ones(n_cells, n_terms);
keepTerm = false(n_terms, 1);
n_kept   = 0;

for k = 1:n_terms
    mask = net.source == terms_all(k);
    tgts = net.target(mask);

    [~, iGene, ~] = intersect(geneNames, tgts, 'stable');
    m_k = numel(iGene);
    if m_k < options.MinGenes
        continue
    end
    keepTerm(k) = true;
    n_kept = n_kept + 1;

    % a: number of top genes that are also in gene set k (per cell)
    a = sum(in_top(iGene, :), 1);   % 1 x n_cells (integer)
    b = n_up - a;                   % top genes not in set
    c_val = m_k - a;                % set genes not in top  (clamp ≥ 0)
    c_val = max(c_val, 0);
    d = n_bg - n_up - c_val;        % neither top nor in set
    d = max(d, 0);

    % Log odds ratio with Haldane-Anscombe correction (eps = 0.5)
    ha = 0.5;
    scores(:, k) = log((a + ha) .* (d + ha) ./ ((b + ha) .* (c_val + ha)))';

    % One-sided Fisher's exact test (enrichment / right tail)
    % p = P(X >= a) where X ~ Hypergeometric(n_bg, m_k, n_up)
    % = 1 - hygecdf(a-1, n_bg, m_k, n_up)
    pv = 1 - hygecdf(a - 1, n_bg, m_k, n_up);
    pvals(:, k) = pv';
end

scores = scores(:, keepTerm);
pvals  = pvals(:, keepTerm);
terms  = terms_all(keepTerm)';

fprintf('SC_ORA: tested %d/%d terms (TopGenes=%d, Background=%d).\n', ...
    n_kept, n_terms, n_up, n_bg);

end
