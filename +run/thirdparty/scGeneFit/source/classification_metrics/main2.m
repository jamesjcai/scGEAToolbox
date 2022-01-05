function benchmarks(data, labels, markers, K)

[d,N]=size(data);
l=labels;

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



baseline=nearest_neighbors_classifier(data, labels, eye(d), test_data, test_labels, K)


P=data(markers,:);
sqz_kmeans=kmeans_classification_error(P, 13, labels)
sqz_knn=nearest_neighbors_classifier(data1, labels1, diag(markers), test_data, test_labels, K)


markers=significant_cluster(G, labels, 'test.png');
aux=zeros(d,1);
aux(markers)=1;
k=5;
nearest_neighbors_classifier(G, labels, diag(aux), test_data, test_labels, k)

P=G(markers,:);

kmeans_classification_error(P, 13, labels)

end