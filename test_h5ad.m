filnm='hca_heart_neuronal_raw.h5ad';
h=h5info(filnm);

% a=h5read(filnm,'/obs/seurat_clusters')
% a=h5read(filnm,'/var/_index')
a=pkg.e_guessh5field(filnm,'/obs/','annotation',true);
a=h5read(filnm,'/obsm/X_umap');
a=h5read(filnm,'/obs/cell_type');