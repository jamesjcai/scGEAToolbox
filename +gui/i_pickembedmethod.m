function [methodtagsel] = i_pickembedmethod(parentfig, allowmulti, dim)

if nargin<3, dim = 0; end
if nargin<2, allowmulti = true; end
if nargin<1, parentfig = []; end

    methodtagsel = [];

    if dim==3
        listitems = {'tSNE 3D', 'UMAP 3D', 'PHATE 3D'};
        methodtag = {'tsne3d', 'umap3d', 'phate3d'};

    elseif dim==2
        listitems = {'tSNE 2D', 'UMAP 2D', 'PHATE 2D'};
        methodtag = {'tsne2d', 'umap2d', 'phate2d'};
    else 
        listitems = {'tSNE 2D', 'tSNE 3D',...
            'UMAP 2D', 'UMAP 3D',...
            'PHATE 2D', 'PHATE 3D'};
        methodtag = {'tsne2d', 'tsne3d', 'umap2d', 'umap3d',...
        'phate2d', 'phate3d', 'metaviz2d', 'metaviz3d'};
    end

    %'MetaViz [PMID:36774377] 2D 🐢',...
    %'MetaViz [PMID:36774377] 3D 🐢'};
    
    
    sce = SingleCellExperiment;
    validmethodtag = fieldnames(sce.struct_cell_embeddings);
    assert(all(ismember(methodtag,validmethodtag)));

    if gui.i_isuifig(parentfig)

        [indx2, tf2] = gui.myListdlg(parentfig, listitems, ...
            'Select embedding methods:',[], allowmulti);
    else
        if allowmulti
            multitag = 'multiple';
        else
            multitag = 'single';
        end
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select embedding methods:'}, ...
                'SelectionMode', multitag, ...
                'ListString', listitems, ...
                'ListSize', [220, 300]);
    end

    if tf2 == 1
        methodtagsel = methodtag(indx2);     
    end
    
end
