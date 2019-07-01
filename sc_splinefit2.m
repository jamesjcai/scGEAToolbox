function [T]=sc_splinefit2(X,Y,genelistx,genelisty,sortid)
if nargin<3, error('X,Y,GENELIST are required.'); end
if nargin<4, genelisty=genelistx; end
if nargin<5, sortid=true; end

[genelist,i,j]=intersect(genelistx,genelisty,'stable');
X=X(i,:);
Y=Y(j,:);

T1=sc_splinefit(X,genelist,false,false);
T1.Properties.VariableNames={'genes','logu1','dropr1','logcv1','d1','pval1','fdr1'};
T2=sc_splinefit(Y,genelist,false,false);
T2.Properties.VariableNames={'genes','logu2','dropr2','logcv2','d2','pval2','fdr2'};

T=join(T1,T2);
T.dd=T.d2-T.d1;
if sortid
    T=sortrows(T,'dd','descend');
end
