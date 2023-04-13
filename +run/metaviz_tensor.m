function [Y]=metaviz_tensor(Sinput,ndim,methodid)


if nargin<3, methodid=1; end
if nargin<2, ndim=2; end

K=length(Sinput);      % K = number of embeddings
n=size(Sinput{1},1);   % n = number of cells
w=zeros(K,n);
%%
D=zeros(n,n,K,'single');
for k=1:K
    d=pdist2(Sinput{k},Sinput{k});
    D(:,:,k)=d./vecnorm(d);
    % (n^2)*(k-1)+1:(n^2)*(k)
end
%%
for x=1:n              % n of cells
    S=1-squareform(pdist(squeeze(D(x,:,:))','cosine'));
    [v,~]=eigs(double(S),1);
    w(:,x)=abs(v);
end

M=zeros(n,n);
for i=1:n             % cell    
    d=zeros(n,1);
    for k=1:K         % type of embedding
        d=d+w(k,i)'.*D(:,i,k);
    end
    M(:,i)=d;
end
M=0.5*(M+M.');

[Y]=pkg.e_embedbyd(M,ndim,methodid);

