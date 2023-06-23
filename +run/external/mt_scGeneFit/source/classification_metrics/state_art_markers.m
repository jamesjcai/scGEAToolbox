load('G2.mat')
load('labels.mat')

[d,N]=size(G);
l=labels;

%split in test and training
rng default
n=ceil(N*0.7);
n_test=N-n;
I=randperm(N,n);
data=G(:,I);
labels1=l(I);
test_data=G;
test_data(:,I)=[];
test_labels=l;
test_labels(I)=[];

markers=distance(G, labels, 'test.png');
aux=zeros(d,1);
aux(markers)=1;
k=5;
nearest_neighbors_classifier(G, labels, diag(aux), test_data, test_labels, k)

P=G(markers,:);

kmeans_classification_error(P, 13, labels)


