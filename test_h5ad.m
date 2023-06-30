filnm='hca_heart_neuronal_raw.h5ad';
h=h5info(filnm);

% a=h5read(filnm,'/obs/seurat_clusters')
% a=h5read(filnm,'/var/_index')
a=h5read(filnm,'/obs/annotation')
a=h5read(filnm,'/obsm/X_umap')
a=h5read(filnm,'/obs/cell_type')