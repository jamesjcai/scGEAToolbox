function [X,g]=e_shadowmatqc(Xmajor,Xminor,gmajor,gminor)

%    [Xmajor,Xminor,gmajor,gminor]=pkg.e_makeshadowmat(sce.X,sce.g);
%    [X1,g1]=pkg.e_shadowmatqc(Xmajor,Xminor,gmajor,gminor);
min_cells_nonzero=20;
% dropout=1;
mtratio=1;
libsize=500;

[X0,g0,keptidxv]=sc_qcfilter(Xmajor,gmajor,libsize,mtratio,min_cells_nonzero);

X1=Xminor;
g1=gminor;
for k=1:length(keptidxv)
    X1=X1(:,keptidxv{k});
end

assert(size(X0,2)==size(X1,2))
X=[X0; X1];
g=[g0; g1];

end



