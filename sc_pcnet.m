function [A]=sc_pcnet(X,ncom)
% [A]=sc_pcnet(X,ncom)
% ncom - number of components used (default=3)
if nargin<2
   ncom=3;
end

% [X]=sc_norm(X);
X=X';
X=zscore(X);

n=size(X,2);
A=1-eye(n);
for k=1:n
    y=X(:,k);
    Xi=X;
    Xi(:,k)=[];
    [~,~,coeff]=svds(Xi,ncom);    
    score=Xi*coeff;
    %[coeff,score]=pca(Xi);
    %coeff=coeff(:,1:ncom);
    score=(score./(vecnorm(score).^2));
    Beta=sum(y.*score);
    A(k,A(k,:)==1)=coeff*Beta';
end

