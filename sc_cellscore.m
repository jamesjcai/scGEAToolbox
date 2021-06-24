function [cs]=sc_cellscore(X,genelist,tgs,nbins,ctrl)
% Create cell scores from a list of genes
%     genes : `list`
%         A list of genes to compute scores from.
%     nbins : `int`
%         Number of expression bins to be used. Default: 25
%     ctrl : `int`
%         Number of control genes in each bin. Default: 100
% Pages 10, 12 and 45 
% https://science.sciencemag.org/content/sci/suppl/2016/04/07/352.6282.189.DC1/Tirosh.SM.pdf
% Line 657
% https://github.com/oscar-franzen/adobo/blob/master/adobo/dr.py

if nargin<5, ctrl=100; end
if nargin<4, nbins=25; end
if nargin<3, [~,tgs]=pkg.i_get_cellcyclegenes; end

X=sc_norm(X);
X=log(X+1.0);
X=zscore(X,0,1);  % https://github.com/satijalab/seurat/issues/1166
X(X>10)=10;   % https://github.com/satijalab/seurat/issues/1166

exprv=pkg.sparse_nanmean(X,2);

rng default
[bin]=discretize(log(10+exprv),nbins);

g=upper(genelist);
[~,idx]=intersect(g,tgs);
i=pkg.i_randsmplbin(idx,bin);

cs=pkg.sparse_nanmean(X(idx,:))-pkg.sparse_nanmean(X(i,:));

end

