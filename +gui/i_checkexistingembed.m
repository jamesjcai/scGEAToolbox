function [vslist] = i_checkexistingembed(sce, ndim)
%see also:             i_pickembedmethod
if nargin<2, ndim=[]; end

vslist = '';
slist = fieldnames(sce.struct_cell_embeddings);
valids = false (length(slist), 1);
for k=1:length(slist)
    sx = sce.struct_cell_embeddings.(slist{k});
    if ~isempty(sx) && size(sx,1) == sce.NumCells
        valids(k) = true;
    end
end
if ~any(valids), return; end
vslist = slist(valids);

if ~isempty(ndim)
   idx = contains(string(vslist), sprintf('%dd',ndim));
   vslist = vslist(idx);
end
