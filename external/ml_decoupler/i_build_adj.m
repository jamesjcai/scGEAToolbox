function [adj, terms] = i_build_adj(geneNames, net, minGenes)
% I_BUILD_ADJ  Build a genes x terms weighted adjacency matrix from a network table.
%
% Shared helper used by sc_mlm, sc_zscore, and other multi-term methods.
% Terms with fewer than minGenes overlapping genes are dropped.
%
% Inputs:
%   geneNames - n_genes x 1 string array (from sce.g)
%   net       - table with columns: source (term), target (gene), weight
%   minGenes  - minimum gene overlap to keep a term
%
% Outputs:
%   adj   - n_genes x n_terms sparse double (weights; 0 = not a member)
%   terms - 1 x n_terms string array of retained term names

if nargin < 3, minGenes = 1; end

geneNames_up = upper(string(geneNames(:)));
terms_all    = unique(net.source, 'stable');
n_terms      = numel(terms_all);
n_genes      = numel(geneNames_up);

net.source = string(net.source);
net.target = upper(string(net.target));

rows = zeros(n_genes * n_terms, 1, 'int32');
cols = zeros(n_genes * n_terms, 1, 'int32');
vals = zeros(n_genes * n_terms, 1);
ptr  = 0;

keepTerm = true(n_terms, 1);
for k = 1:n_terms
    mask = net.source == terms_all(k);
    tgts = net.target(mask);
    wts  = net.weight(mask);

    [~, iGene, iNet] = intersect(geneNames_up, tgts, 'stable');
    if numel(iGene) < minGenes
        keepTerm(k) = false;
        continue
    end
    n_k = numel(iGene);
    rows(ptr+1:ptr+n_k) = int32(iGene);
    cols(ptr+1:ptr+n_k) = int32(k);
    vals(ptr+1:ptr+n_k) = wts(iNet);
    ptr = ptr + n_k;
end

rows = rows(1:ptr);
cols = cols(1:ptr);
vals = vals(1:ptr);

adj_full = sparse(double(rows), double(cols), vals, n_genes, n_terms);
adj   = adj_full(:, keepTerm);
terms = terms_all(keepTerm)';

end
