function [A,coeffv]=sc_pcnetbat(X,ncom)
%Construct network using PC regression
% >>load example_data\pcnet_example.mat
% [X]=log(1+sc_norm(X));     % pcnet input should be LogNormalized
% [A]=sc_pcnet(X,ncom);      % X = expression matrix of genes x cells
% ncom - number of components used (default=3)
% ref: https://rdrr.io/cran/dna/man/PCnet.html
% https://github.com/cran/dna/blob/master/src/rpcnet.c
% https://rdrr.io/cran/dna/f/inst/doc/Introduction.pdf

if nargin<2 || isempty(ncom), ncom=3; end

X=X.';
X=zscore(X);
n=size(X,2);
A=1-eye(n);

D=[];
for k=1:n
    Xi=X;
    Xi(:,k)=[];
    D=cat(3,D,Xi);
end
[~,~,coeffv]=pagesvd(D,"econ");

for k=1:n
    y=X(:,k);
    Xi=D(:,:,k);
    coeff=coeffv(:,1:ncom,k);
    score=Xi*coeff;
    score=(score./(vecnorm(score).^2));
    Beta=sum(y.*score);
    A(k,A(k,:)==1)=coeff*Beta';
end
end
