function [Xmajor, Xminor, gmajor, gminor] = e_makeshadowmat(X, genelist)
%E_MAKESHADOWMAT
%Usage: [Xmajor,Xminor,gmajor,gminor]=pkg.e_makeshadowmat(X,g)
%
% see also: PKG.E_SHADOWMATQC
%
%    [Xmajor,Xminor,gmajor,gminor]=pkg.e_makeshadowmat(sce.X,sce.g);
%    [X1,g1]=pkg.e_shadowmatqc(Xmajor,Xminor,gmajor,gminor);


[T] = sc_genestats(X, genelist);
[~, idx] = maxk(T.Mean, 10);
Xminor = X(idx, :);
gminor = genelist(idx);
Xmajor = X;
Xmajor(idx, :) = [];
gmajor = genelist;
gmajor(idx) = [];
end
