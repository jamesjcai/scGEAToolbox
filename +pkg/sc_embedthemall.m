function [sce] = sc_embedthemall(sce)

if isempty(sce.struct_cell_embeddings.('tsne3d'))
    sce = sce.embedcells('tsne3d',true,true,3);
end

if isempty(sce.struct_cell_embeddings.('tsne2d'))
    sce = sce.embedcells('tsne2d',true,true,2);
end

if isempty(sce.struct_cell_embeddings.('umap3d'))
    sce = sce.embedcells('umap3d',true,true,3);
end

if isempty(sce.struct_cell_embeddings.('umap2d'))
    sce = sce.embedcells('umap2d',true,true,2);
end

if isempty(sce.struct_cell_embeddings.('phate3d'))
       sce = sce.embedcells('phate3d',true,true,3);
end

if isempty(sce.struct_cell_embeddings.('phate2d'))
   sce = sce.embedcells('phate2d',true,true,2);
end

end