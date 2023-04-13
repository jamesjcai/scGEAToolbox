function [Y]=metaviz_memmap(Sinput,ndim,methodid)

if nargin<3, methodid=1; end
if nargin<2, ndim=2; end

K=length(Sinput);      % K = number of embeddings
n=size(Sinput{1},1);   % n = number of cells
w=zeros(K,n);
%%

N=n*n*K;

mmf=tempname;
fileID = fopen(mmf,'w');
for k=1:K
    fwrite(fileID, zeros([n*n,1]),'single');
end
fclose(fileID);
m=memmapfile(mmf,'Format','single','Writable',true);

% D=zeros(n,n,K,'single');
for k=1:K
    d=pdist2(Sinput{k},Sinput{k});
%    D(:,:,k)=d./vecnorm(d);
    m.Data((n^2)*(k-1)+1:(n^2)*(k))=d./vecnorm(d);
end

%isequal(D(:),m.Data)

%%
m.Offset=0;
for x=1:n              % n of cells    
    d=reshape(m.Data(x:n:N),[n K]);
    S=1-squareform(pdist(d','cosine'));

    %S1=1-squareform(pdist(squeeze(D(x,:,:))','cosine'));
    %isequal(S,S1)

    [v,~]=eigs(double(S),1);
    w(:,x)=abs(v);
end

m.Offset=0;
M=zeros(n,n);
for l=1:n             % cell    
    d=zeros(n,1);
    for k=1:K         % type of embedding
        s1=(l-1)*n+1; s2=(n*n)*(k-1); s=s1+s2; t=s+n-1;
        % isequal(s:t,D(:,l,k)')
        d=d+w(k,l)'.*m.Data(s:t);
        %d=d+w(k,l)'.*D(:,l,k);
    end
    M(:,l)=d;
end
M=0.5*(M+M.');

if exist(mmf,'file')==2, delete(mmf); end
[Y]=pkg.e_embedbyd(M,ndim,methodid);

end



