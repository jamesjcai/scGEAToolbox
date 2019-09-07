function [A]=sc_pcnet(X,ncom,fastersvd)
% [A]=sc_pcnet(X,ncom)      % X = expression matrix of genes x cells
% ncom - number of components used (default=3)
% ref: https://rdrr.io/cran/dna/man/PCnet.html
% https://github.com/cran/dna/blob/master/src/rpcnet.c
% https://rdrr.io/cran/dna/f/inst/doc/Introduction.pdf

if nargin<2, ncom=3; end
if nargin<3, fastersvd=false; end
opts.maxit=150;

if fastersvd    
    pw1=fileparts(which(mfilename));
    pth=fullfile(pw1,'thirdparty/faster_svd/lmsvd');
    addpath(pth);
end

%[X]=sc_norm(X);
X=X';
X=zscore(X);

n=size(X,2);
A=1-eye(n);
for k=1:n
    y=X(:,k);
    Xi=X;
    Xi(:,k)=[];
    if fastersvd
        [~,~,coeff]=lmsvd(Xi,ncom,opts);
    else
        [~,~,coeff]=svds(Xi,ncom);
    end
    score=Xi*coeff;
    % [coeff,score]=pca(Xi);
    % coeff=coeff(:,1:ncom);
    score=(score./(vecnorm(score).^2));    
    Beta=sum(y.*score);
    A(k,A(k,:)==1)=coeff*Beta';
end

