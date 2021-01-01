function [ retrievalacc, meanacc] = CalACC1( W, labels )
W = W - diag(diag(W));
U = unique(labels);

if size(labels,1)<size(labels,2)
    labels = labels';
end


for i = 1 : length(U)
    indexx{i} = find(labels==U(i));
    leni(i) = length(indexx{i})-1;
end

[Temp, SortedIndex] = sort(W,2,'descend');

 for i = 1 : length(labels)
     li = labels(i);
     li = find(U==li);
     retrievalacc(i) = length(intersect(indexx{li},SortedIndex(i,1:leni(li))))/leni(li);
 end
     meanacc = mean(retrievalacc);
 end

