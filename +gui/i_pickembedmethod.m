function [methodtagsel] = i_pickembedmethod(parentfig)
if nargin<1, parentfig = []; end

    methodtagsel = [];
    
    listitems = {'tSNE 2D', 'tSNE 3D',...
        'UMAP 2D', 'UMAP 3D',...
        'PHATE 2D', 'PHATE 3D'};
    
    %'MetaViz [PMID:36774377] 2D 🐢',...
    %'MetaViz [PMID:36774377] 3D 🐢'};
    
    methodtag = {'tsne2d', 'tsne3d', 'umap2d', 'umap3d',...
        'phate2d', 'phate3d', 'metaviz2d', 'metaviz3d'};
    
    sce = SingleCellExperiment;
    validmethodtag = fieldnames(sce.struct_cell_embeddings);
    assert(all(ismember(methodtag,validmethodtag)));
    if gui.i_isuifig(parentfig)
        [indx2, tf2] = gui.myListdlg(parentfig, listitems, ...
            'Select embedding methods:');
    else
        [indx2, tf2] = listdlg('PromptString', ...
            {'Select embedding methods:'}, ...
            'SelectionMode', 'multiple', ...
            'ListString', listitems, ...
            'ListSize', [220, 300]);
    end

    if tf2 == 1
        methodtagsel = methodtag(indx2);     
    end
    
end
