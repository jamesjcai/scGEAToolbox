function [markerlist] = sc_pickmarkers(X, genelist, c, ...
    topn, methodid)
if nargin < 5, methodid = 1; end
if nargin < 4, topn = 10; end
assert(isequal(findgroups(c), c));
markerlist = cell(max(c), 1);

switch methodid
    case 1 % Fast method
        % Normalise and log-transform first, the way the sibling marker
        % function does (pkg.e_findallmarkers line 15) and the way the
        % slow method below does with sc_transform. Only this branch --
        % the default, and the one the GUI's Find All Markers and the
        % CLI's --method fast both take -- used to score raw counts.
        %
        % run.ml_PickMarkers ranks by sum_j |mean_i - mean_j|, an
        % absolute difference in count units with no division by
        % expression level or library size, so on raw counts the ranking
        % follows abundance and depth. Measured on a 300-gene fixture
        % with 10 planted markers per cluster and cluster 2 sequenced
        % twice as deep, cluster 2 recovered 0 of its 10 markers -- its
        % top ten filled with ubiquitous high-abundance genes, whose
        % argmax lands in the deepest cluster. After the transform,
        % 10 of 10.
        [idxv] = run.ml_PickMarkers(log1p(sc_norm(X)), genelist, c, topn);
        for k = 1:max(c)
            idx = idxv(1 + (k - 1) * topn:k * topn);
            markerlist{k} = genelist(idx(~isnan(idx)));
        end
    case 2
        num_markers = topn * numel(unique(c));
        markerlist = run.ml_scGeneFit(X, genelist, c, num_markers);
    case 3 % LASSO (slower method)
        for k = 1:max(c)
            fprintf('Processing cell group ... %d of %d\n', k, max(c));
            markerlist{k} = i_pickmarkerslasso(X, genelist, c, k, topn);
        end
    case 4 % Slowest method
        % Transformed once for all groups. I_PICKMARKERS used to call
        % SC_TRANSFORM itself, so the whole matrix was transformed again
        % for every group.
        Xt = sc_transform(X);
        [zRest, zPair] = i_rankzscores(Xt, c);
        for k = 1:max(c)
            a = i_pickmarkers(zRest, zPair, genelist, c, k);
            markerlist{k} = a(1:topn);
        end
    otherwise
        error('sc_pickmarkers:InvalidMethod', 'Unknown methodid %d. Use 1, 2, 3, or 4.', methodid);
end
end


function [markerlist] = i_pickmarkerslasso(X, genelist, idv, id, topn)
idx = idv == id;
y = double(idx);
if issparse(X)
    X = full(X);
end
[B] = lasso(X', y, 'DFmax', topn * 3, 'MaxIter', 1e3);
[~, ix] = min(abs(sum(B > 0) - topn));
b = B(:, ix);
idx = b > 0;
if ~any(idx)
    warning('No marker gene found')
    markerlist = [];
    return;
else
    markerlist = genelist(idx);
    [~, jx] = sort(b(idx), 'descend');
    markerlist = markerlist(jx);
end
end


function [zRest, zPair] = i_rankzscores(X, idv)
% RANKSUM's z for every gene and every comparison I_PICKMARKERS makes.
%   ZREST(:, k)    ranksum(X(:, idv ~= k), X(:, idv == k)), the rest first
%   ZPAIR(:, a, b) ranksum(X(:, idv == a), X(:, idv == b))
% One call covers every group against the rest, and each pair is computed
% once: swapping RANKSUM's two samples only flips the sign of z. This used
% to be a PARFOR over genes calling RANKSUM K^2 times (5 groups x 2000
% genes x 3000 cells: 31.7 s without a parallel pool).
K = max(idv);
[~, ~, z] = pkg.e_ranksumrows(X, idv, "approximate");
zRest = -z;                                  % group first -> rest first
zPair = zeros(size(X, 1), K, K);
for a = 1:K
    for b = a + 1:K
        pair = idv == a | idv == b;
        [~, ~, z] = pkg.e_ranksumrows(X(:, pair), 1 + (idv(pair) == b), "approximate");
        zPair(:, a, b) = z(:, 1);
        zPair(:, b, a) = -z(:, 1);
    end
end
end


function [markerlist, A] = i_pickmarkers(zRest, zPair, genelist, idv, id)
% IDV - cluster ids of cells
% ID  - the id of the cluster, for which marker genes are being identified.
% ZREST, ZPAIR - RANKSUM z-scores from I_RANKZSCORES.
% see also: run.celltypeassignation
% Demo:
% gx=sc_pickmarkers(X,genelist,cluster_id,2);
% run.celltypeassignation(gx)
K = max(idv);
totn = sum(idv ~= id);
A = zeros(size(zRest, 1), K);  % col 1 = all-vs-rest, cols 2..K = per-group
A(:, 1) = zRest(:, id);
col = 1;
for k = 1:K
    if k ~= id
        fprintf('Comparing group #%d with group #%d (out of %d)\n', ...
            id, k, K - 1);
        w = sum(idv == k) ./ totn;
        col = col + 1;
        A(:, col) = w * zPair(:, k, id);     % group k first, as ranksum(x0, x1)
    end
end
A = A(:, 1:col);
[~, idx] = sort(sum(A, 2, 'omitnan'));
markerlist = genelist(idx);
end
