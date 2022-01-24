function [No_cluster,W,idx,eigenvalues,H] = SoptSC_Main(nC,data)
% This function perform the SoptSC algrithm.
%
% Input:
%        nC: Number of cluster.
%        data: A m*n matrix with m rows(genes) and n columns(cells).
%
% Output:
%           W: Cell-to-cell similarity matrix.
%  No_cluster: Number of cluster computed by SoptSC if nC = [];
%             otherwise, No_cluster = nC
%         idx: Cluster label
% eigenvalues: Eigenvalues of graph Laplacian of the consensus matrix
%           H: Non-negative matrix such that W = H*H^T

realdata = data;
realdata = realdata-min(realdata(:));
realdata = realdata./max(realdata(:));

[~,n] = size(realdata);
for i = 1:n
    realdata(:,i) = realdata(:,i)/norm(realdata(:,i));
end

lambda = 0.5;

W = SimilarityM(realdata,lambda,data);

    WB = W;
    n = size(W,1);
    D = diag(WB*ones(n,1));
    Prw = eye(size(W)) - D^(-1/2)*WB*D^(-1/2);
    if n>=1000
        No_eigs = 100;
        all_eigs = real(eigs(Prw,No_eigs,'sm'));
    else
        all_eigs = real(eig(Prw));
    end
    
    ZZ = sort(abs(real(all_eigs)));        
    No_cluster1 = length(find(ZZ<=0.01));
    
%% Determinning the number of clusters
eigenvalues = [];
if isempty(nC)
    [eigenvalues,No_cluster] = Num_cluster(W,No_cluster1);
    nC = No_cluster;
end


flag = 1;
[chuzhiA,~] = nndsvd(W,nC,flag);
params.tol = 10^(-6);
params.Hinit = chuzhiA;
[H,~,~] = symnmf_newton(W, nC, params);
[~,idx] = max(H,[],2);
No_cluster  = nC;
end