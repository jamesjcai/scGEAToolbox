
data_filename='../data/CITEseq.mat';
label_filename='../data/CITEseq-labels.mat';
names_filename='../data/CITEseq_names.mat';

opt.samples_fraction=[0.15];
opt.constraints_neighbors=3;
opt.hinge_scale=7
num_markers=10


aux1=load(data_filename);
aux2=load(label_filename);

tic
[markers,proj_data]=scGeneFit_centers(aux1.G, {aux2.labels}, num_markers, opt);
toc

markers_filename='../output/citeseq_'+string(num_markers)+'markers.csv';
proj_data_filename='../output/citeseq_data_'+string(num_markers)+'markers.csv';
csvwrite(markers_filename,markers);
csvwrite(proj_data_filename,proj_data);

%tsne PLOT
load(names_filename);
plot_filename='../output/citeseq_'+string(num_markers)+'markers.pdf';
rng(100)
tic
P=tsne(proj_data');
toc
plot_data(P, citeseq_names, plot_filename);

%for l=1:max(aux2.zeisel_labels1)
%    plot_filename='../output/zeisel_plot_cluster'+string(l)+'_'+string(num_markers)+'markers.png';
%    rng default 
%    [indices,~]=subcluster_by_label(aux2.zeisel_labels1, l, aux3.zeisel_labels2);
%    plot_data(P, zeisel_names(:,2), plot_filename, indices);
%end
