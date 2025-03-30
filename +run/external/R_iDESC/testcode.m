load ../../../example_data/new_example_sce.mat
% sce2 = run.r_seurat(sce.X, sce.g, './', true);
X = full(sce.X);
save input X -v7.3
writematrix(sce.g, 'g.txt');
t = pkg.makeattributestable(sce);
sequencing_depth = sum(X)';
t = [t table(sequencing_depth)];
writetable(t,'a.txt');

T = readtable('output.csv');

[paramset] = gui.i_degparamset(false, []);
    mindiffpct = paramset{1};
    minabsolfc = paramset{2};
    apvaluecut = paramset{3};
    sortbywhat = paramset{4};
[Tup, Tdn] = pkg.e_processDETable(T,[],[]);