function [A]=sc_pcnetpar(X,ncom)

if nargin<2
   ncom=3;
end

% [X]=sc_norm(X);
X=X';
X=zscore(X);

n=size(X,2);
A=1-eye(n);
B=A(:,1:end-1);

parfor k=1:n
    y=X(:,k);
    Xi=X;
    Xi(:,k)=[];
    [~,~,coeff]=svds(Xi,ncom);
    score=Xi*coeff;
    score=(score./(vecnorm(score).^2));
    Beta=sum(y.*score);
    B(k,:)=coeff*Beta';
end
for k=1:n
   A(k,A(k,:)==1)=B(k,:); 
end
