%% Demonstration of Feature Selection Functions in scGEAToolbox
%% HVG analysis with single data X
%%
cdgea; % set working directory
[X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');   
[X,genelist]=sc_selectg(X,genelist,3,1);
% Normalize data with DESeq method
Xn=sc_norm(X,'type','deseq');
[T]=sc_hvg(Xn,genelist,true,true);

% Highly variable genes (HVGenes), FDR<0.05
HVGenes=T.genes(T.fdr<0.05)

%% Spline-fit feature selection with single data X
%%
[X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');   
[X,genelist]=sc_selectg(X,genelist,3,1);

sortit=true;
[T1]=sc_splinefit(X,genelist,sortit);
% Top 50 featured genes with highest deviation (D) values 
T1.genes(1:50)
dofit=true;
showdata=true;
% Show data points and the spline-fit curve
figure;
sc_scatter3genes(X,genelist,dofit,showdata);
view([36.39 46.25])

%% Analysis of differentially deviated (DD) genes using spline-fit feature selection with data X and Y
%%
% Read and pre-process two data sets, X and Y
[X,genelistx]=sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
[Y,genelisty]=sc_readfile('example_data/GSM3204305_P_N_Expr.csv');
[X,genelistx]=sc_selectg(X,genelistx,3,1);
[Y,genelisty]=sc_selectg(Y,genelisty,3,1);

% Show 3D scatter plot and spline-fit curve for X
figure;
dofit=true;
showdata=true;
subplot(2,1,1)
sc_scatter3genes(X,genelistx,dofit,showdata);
title('Data 1')
view([-6.39 36.70])

% Show 3D scatter plot and spline-fit curve for Y
%figure;
subplot(2,1,2)
sc_scatter3genes(Y,genelisty,dofit,showdata);
title('Data 2')
% view([24.08 32.68])
view([-6.39 36.70])

%% Using function SC_SPLINEFIT2 to fit X and Y separately and obtain DD
% value for each gene.
[T2]=sc_splinefit2(X,Y,genelistx,genelisty,true);

%% Top 50 genes with highest DD value.
T2.genes(1:50)

%% Run GSEAPreranked App with genes ranked with DD
%%
addpath('thirdparty/GSEA');
a=gui.uifig_gseapreranked;
a.load_datatable(T2);

%% Click gene set names in GSEA result table to show GSEA plots


%% The End