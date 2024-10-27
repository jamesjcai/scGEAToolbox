function results = evalCOAL(clusters, n_clusters)
% clusters: each row is a clustering.
[~, m] = size(clusters);
X = zeros(m);
for i = 1:1:m
    X(i, :) = sum(clusters == clusters(:, i));
end
results = runLWEA(X, n_clusters);
end
