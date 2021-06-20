function [T]=sc_genestats(X,genelist)

if isa(X, 'SingleCellExperiment')
    genelist=X.g;
    X=X.X;
else
    if nargin<2, genelist=string(1:size(X,1))'; end
end

dropr=1-sum(X>0,2)./size(X,2);
u=nanmean(X,2);
cv=nanstd(X,[],2)./u;
T=table(genelist,u,cv,dropr);
T.Properties.VariableNames={'Gene','Mean','CV','Dropout_rate'};
% gui.i_viewtable(T);