function [T]=sc_genestat2(X,Y,genelist)

if size(X,1)~=size(Y,1), error('xxx'); end
if nargin<3, genelist=[]; end
[lgu0,dropr0,lgcv0,genelist0]=sc_genestat(X,genelist,false);
[lgu1,dropr1,lgcv1,genelist1]=sc_genestat(Y,genelist,false);
if any(genelist0~=genelist1), error('yyy'); end
T=table(genelist,lgu0,dropr0,lgcv0,lgu1,dropr1,lgcv1);
i=dropr0==1 & dropr1==1;
warning('%d lowly-expressed genes removed.',sum(i));
T(i,:)=[];
%[xyz,i]=sortrows([lgu dropr lgcv],[1 2 3]);

