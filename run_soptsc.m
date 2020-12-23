function [C,W,eigenvalues,H]=run_soptsc(X,varargin)
% run_soptsc - 
%
% REF: SoptSC: Similarity matrix optimization for clustering, lineage, and signaling inference
%
% Similarity matrix-based optimization for single-cell data analysis (SoptSC), 
% in which unsupervised clustering, pseudotemporal ordering, lineage inference, 
% and marker gene identification are inferred via a structured cell-to-cell 
% similarity matrix
%
% Input X: 
% 
% USAGE:
% >> % [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% load('example_data/example10xdata.mat');
% [C]=run_soptsc(X);
% s=sc_tsne(X,3);
% scatter3(s(:,1),s(:,2),s(:,3),20,C,'filled');


if nargin < 1
    error(message('scgea:run_soptsc:TooFewInputs'));
end

p = inputParser;
addRequired(p,'X',@isnumeric);
addOptional(p,'k',[],@isnumeric);
addOptional(p,'donorm',true,@islogical);
parse(p,X,varargin{:});

donorm=p.Results.donorm;
k=p.Results.k;

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/SoptSC');
addpath(pth);
pth=fullfile(pw1,'thirdparty/SoptSC/NNDSVD');
addpath(pth);
pth=fullfile(pw1,'thirdparty/SoptSC/symnmf2');
addpath(pth);

if donorm
    [X]=sc_norm(X,'type','deseq');
%    X=log10(X+1);
end


[coeff, ~, pca_eigvalue] = pca(X');
[~,No_Comps] = max(abs(pca_eigvalue(2:end-1) - pca_eigvalue(3:end)));

aa = max(coeff(:,1:No_Comps+1)');
bb = sort(aa,'descend');
No_features=2000;
if size(X,1) <=1000
    No_sel_genes = size(X,1);
else
    No_sel_genes = min([No_features round(size(X,1))]);
end

gene_selection = aa>=bb(No_sel_genes);
X=X(gene_selection,:);
if isempty(k)
    warning('Number of cluster, k, will be estimated.');
end
[No_cluster,W,C,eigenvalues,H] = SoptSC_Main(k,X);
C=C';

%{
sct=grpstats(s',C,@mean);
figure;
scatter3(s(:,1),s(:,2),s(:,3),20,C,'filled')
hold on
scatter3(sct(:,1),sct(:,2),sct(:,3),180,'r','filled')

figure;
gene_idxv = GC_htmp_DE(X,genelist,cluster_labs,10);

figure; plot_marker(X,{'ACTB','SSR4','PPIB'},genelist,s);
figure; boxplot_marker(X,genelist,{'ACTB','SSR4','PPIB'},C,6);
%}

