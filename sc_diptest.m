function [T]=sc_diptest(X,genelist,sortit,donorm)

% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [X]=sc_norm(X,'type','libsize');
% >> [T]=sc_diptest(X,genelist,true,false);

%ref: https://link.springer.com/chapter/10.1007%2F978-3-319-95933-7_90
% 2.5 Bimodal Genes
% P-values attributed to genes that reject the null hypothesis (unimodal 
% genes) were computed by the dip.test function in R using the default 
% setting. Obtained P-values were adjusted by BH criterion, and genes 
% associated with adjusted P-values less than 0.01 were selected.

if nargin<4, donorm=false; end
if nargin<3, sortit=true; end
if donorm, [X]=sc_norm(X,'type','libsize'); end
u=nanmean(X,2);
cv=nanstd(X,0,2)./u;

n=size(X,1);
%pval=zeros(n,1);
dip=zeros(n,1);
for k=1:n
    %[d,p]=hartigansdipsigniftest(X(k,:),500);
    [d]=hartigansdiptest(X(k,:));
    % pval(k)=p;
    dip(k)=d;
end
% [~,~,~,fdr]=fdr_bh(pval);

% T=table(genelist,u,cv,dip,pval,fdr);
T=table(genelist,u,cv,dip);
T.Properties.VariableNames(1)={'genes'};
if sortit
    T=sortrows(T,'dip','descend');
end

