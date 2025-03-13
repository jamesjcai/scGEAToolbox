function [s] = i_pickembedvalues(sce, ndim, parentfig)
if nargin<2, ndim=[]; end
if nargin<3, parentfig=[]; end

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
if isempty(vslist)
gui.myWarndlg(parentfig, ['No embedding is available. Please embed cells ' ...
        'using an embedding algorithm first.']);
    return; 
end

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, vslist, 'Select an embedding S:');
else
    [indx, tf] = listdlg('PromptString', ...
        {'Select an embedding S:'}, ...
        'SelectionMode', 'single', 'ListString', ...
        vslist, 'ListSize', [220, 300]);
end
if ~tf, return; end
s = sce.struct_cell_embeddings.(vslist{indx});