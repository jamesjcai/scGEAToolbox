function [markerlist,A]=sc_pickmarkers(X,genelist,idv,id)

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
        fprintf('Comparing class %d with class %d (out of %d)\n',...
            id,k,K);
        x0=X(:,idv==k);
        T=i_sc_deg(x0,x1,genelist);
        %a=-log(T.p_val).*sign(T.avg_logFC);
        a=T.z_val;
        A=[A a];        
    end
end
% [~,idx]=sort(sum(A,2));
%A(isnan(A))=0;
%[~,idx]=sort(vecnorm(A,2,2),'descend');  % NaN messed up 
[~,idx]=sort(-vecnorm(A,2,2));  % NaN is ignored
markerlist=genelist(idx);
end


function [T]=i_sc_deg(X,Y,genelist)
ng=size(X,1);
assert(isequal(ng,size(Y,1)));

p_val=ones(ng,1);
avg_logFC=ones(ng,1);
pct_1=ones(ng,1);
pct_2=ones(ng,1);
z_val=ones(ng,1);
for k=1:ng
    x=X(k,:);
    y=Y(k,:);
    [xp,~,xt]=ranksum(x,y);
    p_val(k)=xp;
    z_val(k)=xt.zval;
    avg_logFC(k)=log2(mean(x)./mean(y));
    pct_1(k)=sum(x>0)./length(x);
    pct_2(k)=sum(y>0)./length(y);
end
    p_val_adj = mafdr(p_val,'BHFDR',true);
    sortid=(1:length(genelist))';
    if size(genelist,2)>1 
        gene=genelist';
    else
        gene=genelist;
    end
    T=table(sortid,gene,p_val,avg_logFC,...
        pct_1,pct_2,p_val_adj,z_val);
    % T=sortrows(T,'p_val_adj','ascend');
end
