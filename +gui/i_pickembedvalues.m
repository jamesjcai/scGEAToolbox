function [s] = i_pickembedvalues(sce, ndim)
%see also:             i_pickembedmethod
if nargin<2, ndim=[]; end

s = [];
% slist = fieldnames(sce.struct_cell_embeddings);
% 
% valids = false (length(slist), 1);
% for k=1:length(slist)
%     sx = sce.struct_cell_embeddings.(slist{k});
%     if ~isempty(sx) && size(sx,1) == sce.NumCells
%         valids(k) = true;
%     end
% end
% if ~any(valids), return; end
% vslist = slist(valids);

[vslist] = gui.i_checkexistingembed(sce, ndim);
[indx, tf] = listdlg('PromptString', ...
    {'Select an embedding S:'}, ...
    'SelectionMode', 'single', 'ListString', ...
    vslist, 'ListSize', [175, 130]);
if ~tf, return; end
s = sce.struct_cell_embeddings.(vslist{indx});