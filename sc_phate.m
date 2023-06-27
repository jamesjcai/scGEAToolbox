function s=sc_phate(X,ndim)
%PHATE embedding of cells
% s=sc_phate(X,3);

%see also: SC_TSNT, SC_UMAP
% s_phate=run.mt_PHATE(X,3,true);
% s_umap=run.mt_UMAP(X,3);

if nargin<2, ndim=3; end
% if ~issparse(X), X=sparse(X); end
s=run.mt_PHATE(X,ndim);