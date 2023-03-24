function [T,xyz1]=hvg_splinefit(X,genelist,sortit,removenan)
%SC_SPLINEFIT identify genes with a profile deviated from normal
%
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
% >> [X]=sc_norm(X,'type','libsize');
% >> [T]=sc_splinefit(X,genelist,true,true);

% if nargin<2 || isempty(genelist)
%     genelist=string(1:size(X,1))';
%     genelist=strcat("gene_",genelist);
% end

if nargin<4, removenan=true; end
if nargin<3, sortit=true; end
if nargin<2, genelist=string(1:size(X,1)); end

[lgu,dropr,lgcv,genes]=sc_genestat(X,genelist,sortit,removenan);

% lgu=zscore(lgu);
% dropr=zscore(dropr);
% lgcv=zscore(lgcv);

[~,i]=max(lgcv);
xyz=[lgu dropr lgcv];

[~,j]=sort(pdist2(xyz,xyz(i,:)));
xyz=xyz(j,:)';
lgu=lgu(j);
dropr=dropr(j);
lgcv=lgcv(j);

%xyz=[lgu dropr lgcv]';

s = cumsum([0;sqrt(diff(lgu(:)).^2 + diff(dropr(:)).^2 + diff(lgcv(:)).^2)]);
pp1 = splinefit(s,xyz,15,0.75);
xyz1 = ppval(pp1,s);

D=pdist2(xyz',xyz1');
d=min(D,[],2);
dx=d(d<=quantile(d,0.9));

distFit = fitdist([-dx;dx],'Normal');
pval=normcdf(d,0,distFit.sigma,'upper');
[~,~,~,fdr]=pkg.fdr_bh(pval);

if ~isempty(genes)
   T=table(genes,lgu,dropr,lgcv,d,pval,fdr);   
else
   T=table(lgu,dropr,lgcv,d,pval,fdr);   
end
% 'variablenames',{'Genes','Log10_Mean','Dropout_Rate','Log10_CV','Deviation_3DFeature'});

if sortit, T=sortrows(T,'d','descend'); end

if length(genes)~=length(genelist)
    warning('Output GENES are less than input GENES (some GENES are removed).');
end

