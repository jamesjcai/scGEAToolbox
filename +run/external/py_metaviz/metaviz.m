% pw1=fileparts(mfilename('fullpath'))
% pth=fullfile(pw1,'..','..','thirdparty','PHATE')
%         addpath(pth);


load 1.txt
load 2.txt
load 3.txt
load 4.txt
load 5.txt


Sinput={X1,X2,X3,X4,X5};

%%
K=length(Sinput);      % K = number of embeddings
n=size(Sinput{1},1);   % n = number of cells
w=zeros(K,n);
%%
D=zeros(n,n,K);
for k=1:K
    d=pdist2(Sinput{k},Sinput{k});
    D(:,:,k)=d./vecnorm(d);
    %D(:,:,k)=d;
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
        d=d+w(k,i)'.*squeeze(D(:,i,k));
    end
    M(:,i)=d;
end
M=0.5*(M+M.');
%s=tsne(M,"NumDimensions",2);
load c.txt
%figure; scatter(s(:,1),s(:,2),[],c)

ndim=2;
Y = randmds(M, ndim);
opt = statset('display','iter');
Y = mdscale(M,ndim,'options',opt,'start',Y, ...
    'Criterion','metricstress');

figure; 
scatter(Y(:,1),Y(:,2),[],c);
