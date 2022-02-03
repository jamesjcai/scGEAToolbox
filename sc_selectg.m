function [X,genelist,idx]=sc_selectg(X,genelist,...
                 min_cells_nonzero, isexpressed_cutoff)
%Select genes by expression levels
%
% https://www.nature.com/articles/s41592-018-0254-1
% 11,308 genes with at least 2 cells having more than 
% 4 unique molecular identifier (UMI) reads per cell.

if nargin<4, isexpressed_cutoff=1; end
if nargin<3, min_cells_nonzero=10; end

nc=sum(X>=isexpressed_cutoff,2);
if min_cells_nonzero<1
    idx=nc>=min_cells_nonzero*size(X,2);
else
    idx=nc>=min_cells_nonzero;
end

X=X(idx,:);
if nargout>1
    genelist=genelist(idx);
end


