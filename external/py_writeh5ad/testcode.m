%{
load ../../../example_data/new_example_sce.mat
run.py_writeh5ad(sce, "out.h5ad", './', true);
%}

