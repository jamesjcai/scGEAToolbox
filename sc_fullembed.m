function [sce]=sc_fullembed(sce)


in_fixfield('tsne','tsne3d');
in_fixfield('umap','umap3d');
in_fixfield('phate','phate3d');
in_fixfield('metaviz','metaviz3d');


in_embed('tsne2d',2);
in_embed('tsne3d',3);
in_embed('umap2d',2);
in_embed('umap3d',3);
in_embed('phate2d',2);
in_embed('phate3d',3);
% in_embed('metaviz2d',2);
% in_embed('metaviz3d',3);


    function in_embed(newf, ndim)
        if ~isfield(sce.struct_cell_embeddings,newf) || isempty(sce.struct_cell_embeddings.(newf))
            sce = sce.embedcells(newf, true, true, ndim);
        end   
    end

    function in_fixfield(oldf,newf)
        if ~isfield(sce.struct_cell_embeddings,newf) && isfield(sce.struct_cell_embeddings,oldf)
            if ~isempty(sce.struct_cell_embeddings.(oldf))
                if size(sce.struct_cell_embeddings.(oldf),2) == 3
                    sce.struct_cell_embeddings.(newf) = single(sce.struct_cell_embeddings.(oldf));
                    sce.struct_cell_embeddings = rmfield(sce.struct_cell_embeddings,oldf);
                end
            end
        end
        if isfield(sce.struct_cell_embeddings, oldf)
            if isempty(sce.struct_cell_embeddings.(oldf))
                sce.struct_cell_embeddings = rmfield(sce.struct_cell_embeddings,oldf);
            end
        end
    end

end