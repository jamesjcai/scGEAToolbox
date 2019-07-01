function [lgu,dropr,lgcv,genes]=sc_stat(X,genelist,sortit,removeinf)

if nargin<4, removeinf=true; end
if nargin<3, sortit=true; end
if nargin<2, genelist=[]; end

dropr=1-sum(X>0,2)./size(X,2);
u=nanmean(X,2);
cv=nanstd(X,[],2)./u;
lgu=log(u);
lgcv=log(cv);
genes=genelist;
if sortit
    [xyz,i]=sortrows([lgu dropr lgcv],[1 2 3]);
    lgu=xyz(:,1);
    dropr=xyz(:,2);
    lgcv=xyz(:,3);
    if ~isempty(genelist), genes=genelist(i); end    
end
if removeinf
    i=isnan(lgu) | isinf(lgu) | isnan(lgcv) | isinf(lgcv);
    lgu(i)=[];
    lgcv(i)=[];
    dropr(i)=[];
    if ~isempty(genelist), genes(i)=[]; end
    if length(genes)~=length(genelist)
        warning('Output GENES are less than input GENES (some GENES are removed).');
    end
end