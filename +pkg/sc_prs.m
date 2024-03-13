function [clustered_mat, summary_df] = sc_prs(A, g, Kclusters)

%{
% Binarizing using cutoff and removing unlinked nodes
A(A < Corr_cutoff) = 0;
[i, j, v] = find(A);
A = A(:, any(A, 1));
A = A(any(A, 2), :);

disp(['Correlation matrix created. Shape: ', num2str(size(A))]);

% [i, j, v] = find(data);
edgelist = table(g(j), g(i), v, 'VariableNames', {'gene1', 'gene2', 'Val'});
edgelist = edgelist(edgelist.Val ~= 0, :);

G = graph(edgelist.gene1, edgelist.gene2);

% Calculating giant component - Largest connected component of the graph
disp('Calculating giant component...');
bins = conncomp(G);
[~, idx] = max(histcounts(bins));
Gc = subgraph(G, bins == idx);
graph_gc = Gc;

figure
plot(graph_gc)
disp(['Giant component shape: ', num2str(size(A))]);

% Getting Laplacian matrix
disp('Performing eigen decomposition...');
if ismultigraph(Gc)
    Gc = simplify(Gc);
    figure
    plot(Gc)
end
%}

A = 0.5*(A + A');
A = A - diag(diag(A));
G = graph(A, g);
L = laplacian(G);
n = length(L);


[V, D] = eigs(L, 21, 'smallestabs');
D(1,:) = [];
D(:,1) = [];
V(:,1) = [];
varx = 1 ./ D;
varx(isinf(varx)) = 0;

% Calculating Covariance
cov = V * varx * V';

prs_matrix = cov.^2;
norm_prs_matrix = prs_matrix ./ diag(prs_matrix);



disp('Summarising...');
orf_name = Gc.Nodes.Name;
deg = diag(L);
eff_orig = mean(norm_prs_matrix, 2);
sens_orig = mean(norm_prs_matrix, 1)';

df = table(orf_name, deg, eff_orig, sens_orig);
eigenvector_centr = centrality(Gc, 'eigenvector');
closeness_centr = centrality(Gc, 'closeness');
df.trans = calculate_clustering_coefficient(adjacency(Gc));
df.eigenvec_centr = eigenvector_centr;
df.closeness_centr = closeness_centr;
df.smallest_eigenvec = V(:, 4);


% Clustering
disp('Clustering...');
row_dist = pdist(norm_prs_matrix, 'seuclidean');
% row_dist_sq = squareform(row_dist);
row_linkage = linkage(row_dist, 'ward');

[~, ~, nds_row] = dendrogram(row_linkage,length(norm_prs_matrix));
% leafOrder = optimalleaforder(row_linkage,row_dist);


% row_linkage = [row_linkage, repmat(length(row_dist), size(row_linkage, 1), 1)];
% for i = 1:size(row_linkage, 1)
%     row_linkage(i, 4) = numel(unique([row_linkage(i, 1), row_linkage(i, 2)]));
% end

col_dist = pdist(norm_prs_matrix', 'seuclidean');
% col_dist_sq = squareform(col_dist);
col_linkage = linkage(col_dist, 'ward');

[~, ~, nds_col] = dendrogram(col_linkage,length(norm_prs_matrix));
% nds_col = optimalleaforder(col_linkage,col_dist);

% col_linkage = [col_linkage, repmat(length(col_dist), size(col_linkage, 1), 1)];
% for i = 1:size(col_linkage, 1)
%     col_linkage(i, 4) = numel(unique([col_linkage(i, 1), col_linkage(i, 2)]));
% end

eff_clust = cluster(row_linkage, 'maxclust', Kclusters);
sens_clust = cluster(col_linkage, 'maxclust', Kclusters);

df.eff_clust = eff_clust;
df.sens_clust = sens_clust;

summary_df = sortrows(df, 'eff_orig', 'descend');
clustered_mat = norm_prs_matrix(nds_row, nds_col);

end