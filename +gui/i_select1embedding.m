function [thiss, clabel, sce] = i_select1embedding(sce)
thiss = [];
clabel = '';

listitems = {''};
methodtagv = fieldnames(sce.struct_cell_embeddings);
for k = 1:length(methodtagv)
    methodtag = methodtagv{k};
    if ~isempty(sce.struct_cell_embeddings.(methodtag))
        listitems = [listitems, methodtag];
    end
end
listitems(1) = [];
if isempty(listitems)
    listitems = [listitems, 'Compute tSNE embedding...'];
    %helpdlg('No embedding is available.','');
    %return;
end

[indx2, tf2] = listdlg('PromptString', ...
    {'Select embedding:'}, ...
    'SelectionMode', 'single', 'ListString', listitems);
if tf2 == 1
    clabel = listitems{indx2};
    switch clabel
        case 'Compute tSNE embedding...'
            sce = sce.embedcells('tsne3d', true);
            thiss = sce.struct_cell_embeddings.tsne;
            clabel = 'tsne3d';
        otherwise
            thiss = sce.struct_cell_embeddings.(clabel);
    end
end


end