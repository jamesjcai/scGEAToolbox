function [markerlist]=sc_pickmarkers2(X,genelist,idv,id)
if nargin<4, id=[]; end
% IDV - cluster ids of cells
% ID  - the id of the cluster, for which marker genes are being identified.

if isempty(id)
    K=max(idv);
    markerlist=cell(K,1);
    for k=1:K
        x0=X(:,idv~=k);
        x1=X(:,idv==k);
        [t]=run_mast(x0,x1,genelist);
        markerlist{k}=t.gene;
    end
else
        x0=X(:,idv~=id);
        x1=X(:,idv==id);
        [t]=run_mast(x0,x1,genelist);
        markerlist=t.gene;
end