function [markerlist]=sc_pickmarkers2(X,genelist,idv)

% IDV - cluster ids of cells
% ID  - the id of the cluster, for which marker genes are being identified.

K=max(idv);
markerlist=cell(K,1);
for k=1:K
    x0=X(:,idv~=k);
    x1=X(:,idv==k);
    [t]=run_mast(x0,x1,genelist);
    markerlist{k}=t.gene;
end