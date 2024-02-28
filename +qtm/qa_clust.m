function [c] = qa_clust(s, k)

if nargin < 2 || isempty(k), k = 4; end
[A] = sc_snngraph(s, k);

end