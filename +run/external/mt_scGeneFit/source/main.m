data_filename = '../data/zeisel_data.mat';
label1_filename = '../data/zeisel_labels1.mat';
label2_filename = '../data/zeisel_labels2.mat';
names_filename = '../data/zeisel_names.mat';

%%
opt.samples_fraction = [0.3, 0.9];
opt.constraints_neighbors = 3;
opt.hinge_scale = 0.1;
num_markers = 20;

aux1 = load(data_filename);
aux2 = load(label1_filename);
aux3 = load(label2_filename);

%%
[markers, proj_data] = scGeneFit(aux1.zeisel_data, {aux2.zeisel_labels1, aux3.zeisel_labels2}, num_markers, opt);

markers_filename = '../output/_zeisel_' + string(num_markers) + 'markers.csv';
proj_data_filename = '../output/_zeisel_data_' + string(num_markers) + 'markers.csv';
csvwrite(markers_filename, markers);
csvwrite(proj_data_filename, proj_data);

%tsne PLOTS
load(names_filename);
plot_filename = '../output/zeisel_plot_' + string(num_markers) + 'markers.pdf';
rng(100)
P = tsne(proj_data');
plot_data(P, zeisel_names(:, 1), plot_filename);

for l = 1:max(aux2.zeisel_labels1)
    plot_filename = '../output/zeisel_plot_cluster' + string(l) + '_' + string(num_markers) + 'markers.pdf';
    rng default
    [indices, ~] = subcluster_by_label(aux2.zeisel_labels1, l, aux3.zeisel_labels2);
    plot_data(P, zeisel_names(:, 2), plot_filename, indices);
end
