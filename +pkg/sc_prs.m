function [M, effe_v, sens_v] = sc_prs(A)

arguments
    A {mustBeNumeric, mustBeReal, mustBeEqualSize(A)}
end

%A = 0.5*(A + A');
%A = A - diag(diag(A));
%G = graph(A);
%L = laplacian(G);
[L, Lnorm] = ten.sbe_laplacian_matrix(A);


[V, D] = eigs(L, 21, 'smallestabs');
V(:,1) = [];
d = diag(D);
d(1) = [];                  % or d = d(2:end);
C = V * diag(1 ./ d) * V';  % covariance matrix?
M = C.^2;
M = M ./ diag(M);           % PRS matrix
W = 1 - eye(size(M, 1));
effe_v = mean(M .* W, 2);
sens_v = mean(M .* W, 1)';

% Summarising...
% gname = G.Nodes.Name;
% deg = diag(L);
% T = table(gname, deg, eff_orig, sens_orig);
% eigenvector_centr = centrality(G, 'eigenvector');
% closeness_centr = centrality(G, 'closeness');
% T.trans = calculate_clustering_coefficient(adjacency(G));
% T.eigenvec_centr = eigenvector_centr;
% T.closeness_centr = closeness_centr;
% T.smallest_eigenvec = V(:, 4);

% clustergram

% Clustering ...

%{
row_dist = pdist(M, 'seuclidean');
row_linkage = linkage(row_dist, 'ward');

[~, ~, nds_row] = dendrogram(row_linkage, length(M));
% leafOrder = optimalleaforder(row_linkage,row_dist);
% row_linkage = [row_linkage, repmat(length(row_dist), size(row_linkage, 1), 1)];
% for i = 1:size(row_linkage, 1)
%     row_linkage(i, 4) = numel(unique([row_linkage(i, 1), row_linkage(i, 2)]));
% end

col_dist = pdist(M', 'seuclidean');
% col_dist_sq = squareform(col_dist);
col_linkage = linkage(col_dist, 'ward');

[~, ~, nds_col] = dendrogram(col_linkage,length(M));
% nds_col = optimalleaforder(col_linkage,col_dist);

% col_linkage = [col_linkage, repmat(length(col_dist), size(col_linkage, 1), 1)];
% for i = 1:size(col_linkage, 1)
%     col_linkage(i, 4) = numel(unique([col_linkage(i, 1), col_linkage(i, 2)]));
% end

eff_clust = cluster(row_linkage, 'maxclust', Kclusters);
sens_clust = cluster(col_linkage, 'maxclust', Kclusters);

T.eff_clust = eff_clust;
T.sens_clust = sens_clust;

summary_df = sortrows(T, 'eff_orig', 'descend');
clustered_mat = M(nds_row, nds_col);
%}

end


function mustBeEqualSize(A)
    [m, n] = size(A);
    if ~isequal(m, n)
        eid = 'Size:notEqual';
        msg = 'Input must be a square matrix.';
        error(eid,msg)
    end
end
