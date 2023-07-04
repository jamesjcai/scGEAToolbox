function [hFig]=i_violinmatrix(X,g,c,~,tgene,~)

if nargin<6, uselog=false; end
[yes]=ismember(tgene,g);
if ~any(yes), return; end
z=length(tgene)-sum(yes);
if z>0
    fprintf('%d gene(s) not in the list are excluded.\n',z); 
end
tgene=tgene(yes);

M=length(tgene);
% N=length(cL);

hFig=figure;
ct=1;
for k=1:M
    a1=X(g==tgene(k),:);
    subplot(M,1,ct)        
    pkg.violinplot(a1.',c,'showdata',false);
    ylabel(tgene(k));
    ct=ct+1;
end

