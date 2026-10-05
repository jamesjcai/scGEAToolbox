%[text] # Demo 2 - Feature Selection Functions
%[text] ## HVG analysis with single data X
cdgea; % set working directory
[X,genelist]=sc_readfile(pkg.i_exampledata('GSM3044891_GeneExp.UMIs.10X1.txt'));
[X,genelist]=sc_selectg(X,genelist,3,1);
% Normalize data with DESeq method
Xn=sc_norm(X,'type','deseq');
[T]=sc_hvg(Xn,genelist,true,true);

% Highly variable genes (HVGenes), FDR<0.05
HVGenes=T.genes(T.fdr<0.05);
disp(HVGenes(1:10))
%%
%[text] ## Spline-fit feature selection with single data X
[X,genelist]=sc_readfile(pkg.i_exampledata('GSM3044891_GeneExp.UMIs.10X1.txt'));
[X,genelist]=sc_selectg(X,genelist,3,1);

sortit=true;
[T1]=sc_splinefit(X,genelist,sortit);
% Top 10 featured genes with highest deviation (D) values 
T1.genes(1:10)
dofit=true;
showdata=true;
% Show data points and the spline-fit curve
gui.i_hvgcurveplot(X,genelist,dofit,showdata,[],"splinefit");
%view([36.39 46.25])
%%
%[text] ## Analysis of differentially deviated (DD) genes using spline-fit feature selection with data X and Y
%[text] Read and pre-process two data sets, X and Y
[X,genelistx]=sc_readfile(pkg.i_exampledata('GSM3204304_P_P_Expr.csv'));
[Y,genelisty]=sc_readfile(pkg.i_exampledata('GSM3204305_P_N_Expr.csv'));
[X,genelistx]=sc_selectg(X,genelistx,3,1);
[Y,genelisty]=sc_selectg(Y,genelisty,3,1);

% Show 3D scatter plot and spline-fit curve for X
% figure;
dofit=true;
showdata=true;
%subplot(2,1,1)
gui.i_hvgcurveplot(X,genelistx,dofit,showdata,[],"splinefit");
title('Data 1')
%view([-6.39 36.70])

% Show 3D scatter plot and spline-fit curve for Y
%figure;
%subplot(2,1,2)
gui.i_hvgcurveplot(Y,genelisty,dofit,showdata,[],"splinefit");
title('Data 2')
% view([24.08 32.68])
% view([-6.39 36.70])
%%
%[text] ## Using function SC\_SPLINEFIT2 to fit X and Y separately and obtain DD value of each gene
[T2]=sc_splinefit2(X,Y,genelistx,genelisty,true);
%%
%[text] ## Top 10 genes with highest DD value.
T2.genes(1:10)
%%
%[text] ## The End

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
