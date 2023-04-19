function [T,sx,sy,sz,d]=sc_splinefit2(X,Y,genelistx,genelisty,sortid)
%Compare genes' expression profiles with spline fit regression
if nargin<3, error('X,Y,GENELIST are required.'); end
if nargin<4, genelisty=genelistx; end
if nargin<5, sortid=true; end

[genelist,i,j]=intersect(genelistx,genelisty,'stable');
X=X(i,:);
Y=Y(j,:);

[T1,sx]=sc_splinefit(X,genelist);
T1.Properties.VariableNames={'genes','logu1','logcv1','dropr1','d1','pval1','fdr1'};
[T2,sy]=sc_splinefit(Y,genelist);
T2.Properties.VariableNames={'genes','logu2','logcv2','dropr2','d2','pval2','fdr2'};

T=join(T1,T2,'Keys','genes');
T.dd=T.d2-T.d1;
if sortid
    T=sortrows(T,'dd','descend');
end

[d,sz] = procrustes(sx,sy);


