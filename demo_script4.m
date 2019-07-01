%% Demonstration of Clustering Functions in scGEApp
%% Load example data
%%
cdgea; % set working directory
% load('example_data/example10xdata2.mat','X','genelist');
[X,genelistx]=sc_readfile('example_data/GSM3204304_P_P_Expr_999cells.csv');
[Y,genelisty]=sc_readfile('example_data/GSM3204305_P_N_Expr_999cells.csv');
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
C=sc_cluster([X Y],'type','simlr');

%% Plot the clustering result
%%
s=sc_tsne([X Y],2,false,true,false);   % s=sc_tsne(X,ndim,plotit,donorm,dolog1p);
%%
figure;
scatter(s(:,1),s(:,2),20,C,'filled')

figure;
scatter(s(:,1),s(:,2),20,cellidx,'filled')

%% The End