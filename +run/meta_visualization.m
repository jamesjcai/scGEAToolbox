function [Y]=meta_visualization(X,ndim)
if nargin<2, ndim=2; end
S=[];
pw1=fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)    
    pth=fullfile(pw1,'thirdparty','PHATE'); % for calling randmds.m
    addpath(pth);
    pth1=fullfile(pw1,'thirdparty','umapFileExchange');
    pth3=fullfile(pw1,'thirdparty','umapFileExchange','umap.jar');
    addpath(pth1);
    javaaddpath(pth3);
end

nstep=6+1;

fw=gui.gui_waitbar_adv;
Xn=log(1+sc_norm(X))';
data = svdpca(Xn, 50, 'random');

gui.gui_waitbar_adv(fw,1/nstep,'Meta Visualization - PCA...');
[~,S{1}]=pca(data,NumComponents=ndim);

gui.gui_waitbar_adv(fw,2/nstep,'Meta Visualization - MDS...');
D=squareform(pdist(data));
S{end+1}=pkg.e_embedbyd(D,ndim,2);

gui.gui_waitbar_adv(fw,3/nstep,'Meta Visualization - TSNE 1/3...');
S{end+1}=tsne(data,Perplexity=30,NumDimensions=ndim);

gui.gui_waitbar_adv(fw,3/nstep,'Meta Visualization - TSNE 2/3...');
S{end+1}=tsne(data,Perplexity=15,NumDimensions=ndim);

gui.gui_waitbar_adv(fw,3/nstep,'Meta Visualization - TSNE 3/3...');
S{end+1}=tsne(data,Perplexity=50,NumDimensions=ndim);

gui.gui_waitbar_adv(fw,4/nstep,'Meta Visualization - UMAP 1/3...');
S{end+1}=run_umap_main(data,'n_components',ndim, ...
    'n_neighbors',15,'verbose','none');

gui.gui_waitbar_adv(fw,4/nstep,'Meta Visualization - UMAP 2/3...');
S{end+1}=run_umap_main(data,'n_components',ndim, ...
    'n_neighbors',30,'verbose','none');

gui.gui_waitbar_adv(fw,4/nstep,'Meta Visualization - UMAP 3/3...');
S{end+1}=run_umap_main(data,'n_components',ndim, ...
    'n_neighbors',50,'verbose','none');

gui.gui_waitbar_adv(fw,5/nstep,'Meta Visualization - PHATE 1/3...');
S{end+1}=phate(sqrt(Xn), 't', 20, 'ndim', ndim, 'k', 5);
gui.gui_waitbar_adv(fw,5/nstep,'Meta Visualization - PHATE 2/3...');
S{end+1}=phate(sqrt(Xn), 't', 20, 'ndim', ndim, 'k', 15, 'pot_method', 'sqrt');
gui.gui_waitbar_adv(fw,5/nstep,'Meta Visualization - PHATE 3/3...');
S{end+1}=phate(sqrt(Xn), 't', 20, 'ndim', ndim, 'k', 30);

gui.gui_waitbar_adv(fw,6/nstep,'Meta Visualization - METAVIZ');
[Y]=run.metaviz(S,ndim);
gui.gui_waitbar_adv(fw);
end

%figure; scatter(Y(:,1),Y(:,2));

%{
disp('sTSNE2')
sTSNE2=tsne(Xn,Perplexity=15,NumDimensions=3);
disp('sTSNE3')
sTSNE3=tsne(Xn,Perplexity=50,NumDimensions=3);
disp('sTSNE4')
sTSNE4=tsne(data,Perplexity=30,NumDimensions=3);
disp('sTSNE5')
sTSNE5=tsne(data,Perplexity=15,NumDimensions=3);
disp('sTSNE6')
sTSNE6=tsne(data,Perplexity=50,NumDimensions=3);


disp('sUMAP2')
sUMAP2=run_umap_main(Xn,'n_components',3,'n_neighbors',30);
disp('sUMAP3')
sUMAP3=run_umap_main(Xn,'n_components',3,'n_neighbors',50);
disp('sUMAP4')
sUMAP4=run_umap_main(data,'n_components',3,'n_neighbors',15);
disp('sUMAP5')
sUMAP5=run_umap_main(data,'n_components',3,'n_neighbors',30);
disp('sUMAP6')
sUMAP6=run_umap_main(data,'n_components',3,'n_neighbors',50);

disp('sPHATE2')
sPHATE2 = phate(sqrt(Xn), 't', 20, 'ndim', 3, 'k', 30);
disp('sPHATE3')
sPHATE3 = phate(sqrt(Xn), 't', 20, 'ndim', 3, 'k', 50);
disp('sPHATE4')
sPHATE4 = phate(data, 't', 20, 'ndim', 3, 'k', 5);
disp('sPHATE5')
sPHATE5 = phate(data, 't', 20, 'ndim', 3, 'k', 30);
disp('sPHATE6')
sPHATE6 = phate(data, 't', 20, 'ndim', 3, 'k', 50);

%}

% PCA: the fast SVD function svds from R package rARPACK with embedding dimension k=2.
% MDS: the basic R function cmdscale with embedding dimension k=2.
% Sammon: the R function sammon from R package MASS with embedding dimension k=2.
% LLE: the R function lle from R package lle with parameters m=2, k=20, reg=2.
% HLLE: the R function embed from R package dimRed with parameters method="HLLE",knn=20, ndim=2.
% Isomap: the R function embed from R package dimRed with parameters method="Isomap",knn=20, ndim=2.
% kPCA1&2: the R function embed from R package dimRed with parameters method="kPCA",kpar=list(sigma=width), ndim=2, where we set width=0.01 for kPCA1 and width=0.001 for kPCA2.
% LEIM: the R function embed from R package dimRed with parameters ndim=2 and method ="LaplacianEigenmaps".
% UMAP1&2: the R function umap from R package uwot with parameters n neighbors=n,n components=2, where we set n=30 for UMAP1 and width=50 for UMAP2.
% tSNE1&2: the R function embed from R package dimRed with parameters method="tSNE",perplexity=n, ndim=2, where we set n=10 for tSNE1 and n=50 for tSNE2.
% PHATE1&2: the R function phate from R package phateR with parameters knn=n, ndim=2,where we set n=30 for PHATE1 and n=50 for PHATE2.
