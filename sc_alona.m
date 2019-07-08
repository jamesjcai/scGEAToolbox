function [cluster_labs,s]=sc_alona(X,genelist,donorm)
% sc_alona - 
%
% This file contains clustering methods used by alona.
%  In general, it flows like this:
%     1. identify highly variable genes (HVG), retrieve N genes
%     2. perform PCA on the HVG, retrieve N components
%     3. adjust PCAs by weight
%     4. compute KNN
%     5. compute SNN from KNN, prune SNN graph
%     6. identify communities with leiden algo
%     7. run t-SNE on the PCAs
%  How to use alona:
%  https://github.com/oscar-franzen/alona/
%  https://alona.panglaodb.se/
%  https://github.com/oscar-franzen/alona/blob/master/alona/clustering.py
 
% USAGE:
% >> % [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% load('example_data/example10xdata.mat');
% [C,s]=sc_alona(X,[],true);
% figure;
% scatter(s(:,1),s(:,2),20,C,'filled')

if nargin<3
    donorm=true;
end
[~,Xsorted]=sc_hvg(X,genelist,true,false,donorm);

X=Xsorted(1:500,:);
% genelist=genelistsorted(1:500);
[~,score]=pca(X');


tic; [idx]=knnsearch(score,score,'K',11,'distance','cosine'); toc;

% http://mlwiki.org/index.php/SNN_Clustering


%         Cluster the SNN graph using the Leiden algorithm.
%         https://github.com/vtraag/leidenalg

s=tsne(X','NumDimensions',ndim,'InitialY',score(:,1:ndim));





