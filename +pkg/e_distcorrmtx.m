function [R]=e_distcorrmtx(X)

n=size(X,1);  % number of genes
R=zeros(n);
c=0;
for k=1:n
    for l=k+1:n
        c=c+1;
        R(k,l)=pkg.distcorr(X(k,:)',X(l,:)');
        %R(k,l)=pkg.bcdistcorr(X(k,:)',X(l,:)');
        if c>20, break; end
    end
        if c>20, break; end

end
