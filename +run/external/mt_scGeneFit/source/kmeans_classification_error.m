function min_error= kmeans_classification_error(P, k, labels)
[d,N]=size(P);

%min over 10 experiments
min_error=100;
for tt=1:10
idx=kmeans(P',k);
error=0;

for i=1:N
    for j=i+1:N
        if and(labels(i)==labels(j), idx(i)~=idx(j))
            error=error+1;
        elseif and(labels(i)~=labels(j), idx(i)==idx(j))
            error=error+1;
        end
    end
end

error=error/(N*(N-1)/2)*100;
if min_error>error
    min_error=error;
end
end