
%Load data
data_filename='../data/zeisel_data.mat';
label1_filename='../data/zeisel_labels1.mat';
label2_filename='../data/zeisel_labels2.mat';
names_filename='../data/zeisel_names.mat';

opt.samples_fraction=[0.1, 0.2];
opt.hinge_scale=1;
num_markers=30;

aux1=load(data_filename);
aux2=load(label1_filename);
aux3=load(label2_filename);

%Compute markers
tic
[markers,proj_data]=scGeneFit_centers(aux1.zeisel_data/100, {aux2.zeisel_labels1, aux3.zeisel_labels2}, num_markers, opt);
toc

%Write the results
markers_filename='../output/5new_zeisel_'+string(num_markers)+'markers.csv';
proj_data_filename='../output/5new_zeisel_data_'+string(num_markers)+'markers.csv';
csvwrite(markers_filename,markers);
csvwrite(proj_data_filename,proj_data);

%tsne PLOTS
load(names_filename);
plot_filename='../output/5new_zeisel_'+string(num_markers)+'markers.pdf';
rng default
P=tsne(proj_data');
plot_data(P, zeisel_names(:,1), plot_filename);
%%
for l=1:max(aux2.zeisel_labels1)
    plot_filename='../output/5new_zeisel_plot_cluster'+string(l)+'_'+string(num_markers)+'markers.pdf';
    rng default 
    [indices,~]=subcluster_by_label(aux2.zeisel_labels1, l, aux3.zeisel_labels2);
    plot_data(P, zeisel_names(:,2), plot_filename, indices);
    
    Q=tsne(proj_data(:,indices)');
    plot_filename='../output/5new_zeisel_plot_subcluster'+string(l)+'_'+string(num_markers)+'markers.pdf';
    rng default
    plot_data(Q, zeisel_names(indices,2), plot_filename);
    
end

%Benchmarks 6,5,2,7
for T=6:6
    filename='../output/significant'+string(T)+'.pdf';
    CLUSTER=T;
[indices,labels]=subcluster_by_label(aux2.zeisel_labels1, T, aux3.zeisel_labels2);
benchmarks_and_plot(aux1.zeisel_data(:,indices)/100, labels, markers, [3,5,15], zeisel_names(indices,2), filename);
end