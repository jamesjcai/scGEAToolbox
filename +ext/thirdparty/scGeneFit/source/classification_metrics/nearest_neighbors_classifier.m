function e= nearest_neighbors_classifier(data, labels, P, test_data, test_labels, k)
%P is dxd
%data is dxn
%test is dxm
%k is list of integers

squeezed= P*data;
squeezed_test= P*test_data;

K=max(k);
I= knnsearch(squeezed', squeezed_test', 'K', K);

%compute errors
e=zeros(size(k));

for i=1:size(test_data,2)
aux=labels(I(i,:));

for t=1:size(k,1)
    tt=k(t);
    if mode(aux(1:tt))~= test_labels(i)
        e(t)=e(t)+1;
    end
end
end

e=e*100/size(test_data,2);
end