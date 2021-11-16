function [F]=knk4_buildPerturbationLandscape_quick(A,genelist,oldF,a,b)

if nargin<5, b=size(A,1); end
if nargin<4, a=1; end
if nargin<3, oldF=[]; end
import ten.*
n=length(genelist);
assert(n==size(A,1));
if ~isempty(oldF)
    F=oldF;
else
    F=zeros(n,n);
end

    ndim=30;
    mu=0.9;
    W1=A+1;
    W2=W1;
    assert(n==size(W1,2));
    W12=eye(size(W1,2),size(W2,2));
    mu = mu*(sum(W1(:))+sum(W2(:)))/(2*sum(W12(:)));
    W = [W1 mu*W12; mu*W12' W2];

    a=min([a n]); b=min([b n]);
    for k=a:b
        fprintf('%s ... gene %d of %d\n',genelist(k),k,n);
        if all(F(:,k)==0)
            tic;
            if k>1
                W(n+k-1,n+1:end)=W(k-1,1:n);
            end
            W(n+k,n+1:end)=1;
            D=sum(abs(W));
            L=diag(D)-W;
            [V,D] = eigs(L,ndim*2,'smallestreal');
            d=diag(D);
            [d,ind] = sort(d);
            V=V(:,ind);
            V=V(:,d>=1e-8);
            V=V(:,1:ndim);
            aln0=V(1:n,:);
            aln1=V(n+1:end,:);
            drdist=vecnorm(aln0-aln1,2,2).^2;
            F(:,k)=drdist./norm(drdist);
            toc;
        end
    end
end
