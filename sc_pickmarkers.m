function [markerlist]=sc_pickmarkers(X,genelist,idv,id)

% IDV - cluster ids of cells
% ID  - the id of the cluster, for which marker genes are being identified.
% see also: run_celltypeassignation
% Demo:
%gx=sc_pickmarkers(X,genelist,cluster_id,2);
%run_celltypeassignation(gx)
K=max(idv);
x1=X(:,idv==id);
A=[];
for k=1:K    
    if k~=id
        x0=X(:,idv==k);
        T=sc_deg(x0,x1,genelist,1);
        %T=sortrows(T,'sortid','ascend');
        a=-log(T.p_val).*sign(T.avg_logFC);
        A=[A a];
    end
end
[~,idx]=sort(sum(A,2));
markerlist=genelist(idx);

