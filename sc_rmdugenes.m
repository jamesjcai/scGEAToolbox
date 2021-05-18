function [X,genelist]=sc_rmdugenes(X,genelist)
% remove genes with duplicate name

[~, w] = unique( genelist, 'stable' );
duplicate_indices = setdiff( 1:numel(genelist), w );
if isempty(duplicate_indices)
    return;
end

for k=1:length(duplicate_indices)
    idx=find(genelist==genelist(duplicate_indices(k)));
    X(idx(1),:)=sum(X(idx,:),1);    
end
X(duplicate_indices,:)=[];
genelist(duplicate_indices)=[];
end

