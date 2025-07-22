load ../../../example_data/new_example_sce.mat
sce2 = run.r_seurat(sce.X, sce.g, './', true);
