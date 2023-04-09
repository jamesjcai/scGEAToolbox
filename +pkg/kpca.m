function [coeff,score,latent]=kpca(X,k,sgm,isds)

% how to use coeff of kpca:
% 	[coeff] = pkg.kpca(XTaining);
%   Ktest = constructKernel(XTest)
%   Y = Ktest*coeff;
%   Y is the embedding result of XTest.
%
% (c) 2023 James Cai jcai@tamu.edu

if nargin<2, k=2; end
if nargin<3, sgm=40; end
if nargin<4, isds=false; end    % if true, then input X is a distance matrix

n=size(X,1);
if isds
    D=X;
else
    try
        d=dot(X,X,2);
        D=d+d'-2*(X*X');
    catch
        D=pdist2(X,X).^2;
    end
end
K=exp(-D/(2*sgm^2));
H=eye(n)-ones(n,n)/n;
Kc=H*K*H;
[coeff,latent]=eigs(Kc,k);   
score=Kc*coeff;
