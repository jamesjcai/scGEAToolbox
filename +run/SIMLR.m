function [C,s]=SIMLR(X,k,donorm)
% RUN_SIMLR - 
%
% REF: https://www.nature.com/articles/nmeth.4207
% Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning
%
% Input X: 
% 
% USAGE:
% >> % [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% load('example_data/example10xdata.mat');
% [C,s]=run_simlr(X,[],true);
% figure;
% scatter(s(:,1),s(:,2),20,C,'filled')

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/SIMLR');
addpath(pth);
pth=fullfile(pw1,'thirdparty/SIMLR/src');
addpath(pth);

if nargin<2
    k=[];
end

if nargin<3
    donorm=true;
end

if donorm
    [X]=sc_norm(X,'type','deseq');
    X=log2(X+1);
end

if isempty(k)
    [K1, K2] = Estimate_Number_of_Clusters_SIMLR(X',2:10);
    [~,i]=min(K2);
    k=i+1;
end
[C,S,F,s,alpha] = SIMLR_ori(X',k,10,0,0);
C=C';

% [C, S, F, s,alpha] = SIMLR_pearson(X',k,10,0,0);
%    figure;
%    gscatter(ydata(:,1),ydata(:,2),y);
% figure;
% scatter(s(:,1),s(:,2),20,C,'filled')

