cdgea;
load example_data\gettingstarted_test.mat X g s
sc_scatter(X,g,s);
pause(1.5)
load example_data\gettingstarted_test_sce sce
sc_scatter_sce(sce);
pause(1.5)
figure;
run.sc_marker(X,g,grp2idx(sce.c),'plotit',true);