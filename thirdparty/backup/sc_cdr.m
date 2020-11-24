function [r]=sc_cdr(X,cutoff)
%SC_CDR - calcualtes Cellular Detection Rate (CDR) of cells
% CDR is the proportion of genes detected in each cell.
% https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5

if nargin<2, cutoff=0; end
% n=size(X,1);   % number of genes
% r=sum(X>cutoff)./n;
r=mean(X>cutoff);

% Usage:
% [~,s]=pca(X');
% scatter(sc_cdr(X),s(:,1))
% xlabel('Cellular Detection Rate')
% ylabel('PC_1')
