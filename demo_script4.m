%% Demo 4 - Clustering Functions
%% Load example data
%%
cdgea; % set working directory
% load('example_data/example10xdata2.mat','X','genelist');
[X,genelistx]=sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
[Y,genelisty]=sc_readfile('example_data/GSM3204305_P_N_Expr.csv');
[X,genelistx]=sc_selectg(X,genelistx,3,1);
[Y,genelisty]=sc_selectg(Y,genelisty,3,1);

%% Intersection of common genes in X, Y and Z
%%
[genelist]=intersect(genelistx,genelisty,'stable');
% Remove genes encoded in the mitochondrial genome
i=startsWith(genelist,'MT-');
genelist(i)=[];
[~,i1]=ismember(genelist,genelistx);
[~,i2]=ismember(genelist,genelisty);
X=X(i1,:); genelistx=genelist;
Y=Y(i2,:); genelisty=genelist;

%% Label cells
%
cellidx=[1*ones(size(X,2),1); 2*ones(size(Y,2),1)];
%% Cluster cells using SIMLR
%%
C=sc_cluster_x([X Y],5,'type','simlr');

% [C,s]=run.SIMLR(X,[]);  % auto-determine the number of clusters.
%% Plot the clustering result
%%
s=sc_tsne([X Y],2,false,true,false);   % s=sc_tsne(X,ndim,plotit,donorm,dolog1p);
%%
figure;
scatter(s(:,1),s(:,2),20,C,'filled')

figure;
scatter(s(:,1),s(:,2),20,cellidx,'filled')

%% Using SC3 example data yan.csv

[X,genelist]=sc_readtsvfile('example_data/yan.csv');
t=readtable('example_data\yan_celltype.txt');
celltypelist=string(t.cell_type1);
rng(235); showlegend=true;

s=sc_tsne(X,2);
c1=run.SC3(X,6);
c2=run.SIMLR(X,6);
c3=run.SoptSC(X,'k',6);

% Result of SC3/R pacakge
load example_data/sc3_results.txt
c0=sc3_results;

%% Compare clustering results between SC3/R vs SC3, SIMILR and SoptSC
Cal_NMI(c0,c1)
Cal_NMI(c0,c2)
Cal_NMI(c0,c3)


fh=figure; 
subplot(2,2,1)
gscatter(s(:,1),s(:,2),celltypelist)
if showlegend, legend('Location','northwest'); else, legend off; end
title('Cell type (''Biological Turth'')')

subplot(2,2,2)
gscatter(s(:,1),s(:,2),c1)
if showlegend, legend('Location','northwest'); else, legend off; end
title('SC3')

subplot(2,2,3)
gscatter(s(:,1),s(:,2),c2)
if showlegend, legend('Location','northwest'); else, legend off; end
title('SIMLR')

subplot(2,2,4)
gscatter(s(:,1),s(:,2),c3)
if showlegend, legend('Location','northwest'); else, legend off; end
title('SoptSC')
fh.Position=[fh.Position(1) fh.Position(2)-100 fh.Position(3)+100 fh.Position(4)+100];


%% The End