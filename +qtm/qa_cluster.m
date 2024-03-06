function [c] = qa_cluster(s, k)

if nargin < 2 || isempty(k), k = 4; end
[A] = sc_snngraph(s, k);
squareform(pdist(x','jaccard'));

% https://doi.org/10.1093/bib/bbad377
end