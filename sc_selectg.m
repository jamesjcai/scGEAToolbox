function [X,genelist,idx]=sc_selectg(X,genelist,min_count,min_cells)
%
% https://www.nature.com/articles/s41592-018-0254-1
% 11,308 genes with at least 2 cells having more than 4 unique molecular identifier (UMI) reads per cell. 
if nargin<4, min_cells=10; end
if nargin<3, min_count=1; end

nc=sum(X>min_count,2);
if min_cells<1
    idx=nc>=min_cells*size(X,2);
else
    idx=nc>=min_cells;
end

X=X(idx,:);
if nargout>1
    genelist=genelist(idx);
end


