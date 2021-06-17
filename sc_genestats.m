function [T]=sc_genestats(X,genelist)

if nargin<2, genelist=string(1:size(X,1))'; end

dropr=1-sum(X>0,2)./size(X,2);
u=nanmean(X,2);
cv=nanstd(X,[],2)./u;
T=table(genelist,u,cv,dropr);
T.Properties.VariableNames={'Gene','Mean','CV','Dropout_rate'};
% fig = uifigure;
% uit = uitable(fig,'Data',T);
