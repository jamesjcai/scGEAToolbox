function benchmarks(data, labels, markers, K)

[d,N]=size(data);
l=labels;
N_clusters=size(unique(labels),1);

%split in 30% test and 70% training
rng default
n=ceil(N*0.7);
n_test=N-n;
I=randperm(N,n);
data1=data(:,I);
labels1=l(I);
test_data=data;
test_data(:,I)=[];
test_labels=l;
test_labels(I)=[];

baseline_kmeans=kmeans_classification_error(data, N_clusters, labels)
K
baseline_knn=nearest_neighbors_classifier(data, labels, eye(d), test_data, test_labels, K)


P=data(markers,:);
sqz_kmeans=kmeans_classification_error(P, N_clusters, labels)
sqz_knn=nearest_neighbors_classifier(data1, labels1, diag(markers), test_data, test_labels, K)


markers_sign=significant_cluster(data, labels);
aux=zeros(d,1);
aux(markers_sign)=1;
significant_knn=nearest_neighbors_classifier(data1, labels1, diag(aux), test_data, test_labels, K)

P=data(markers_sign,:);
significant_kmeans=kmeans_classification_error(P, N_clusters, labels)

markers_rand=randperm(d,40);
aux=zeros(d,1);
aux(markers_rand)=1;
rand40_knn=nearest_neighbors_classifier(data1, labels1, diag(aux), test_data, test_labels, K)

P=data(markers_rand,:);
rand40_kmeans=kmeans_classification_error(P, N_clusters, labels)


markers_rand=randperm(d,7);
aux=zeros(d,1);
aux(markers_rand)=1;
rand7_knn=nearest_neighbors_classifier(data1, labels1, diag(aux), test_data, test_labels, K)

P=data(markers_rand,:);
rand7_kmeans=kmeans_classification_error(P, N_clusters, labels)
end