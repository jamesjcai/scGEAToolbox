function s=sc_umap(X,ndim)
%UMAP embedding of cells
% s=sc_umap(X,3);

%see also: SC_TSNT, SC_PHATE
% s_phate=run.PHATE(X,3,true);
% s_umap=run.UMAP(X,3);

if nargin<2, ndim=3; end
% if ~issparse(X), X=sparse(X); end
s=run.UMAP(X,ndim);