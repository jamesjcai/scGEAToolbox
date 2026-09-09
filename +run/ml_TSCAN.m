function [t] = ml_TSCAN(X, varargin)
%ML_TSCAN  Pseudotime ordering by TSCAN (Ji & Ji, NAR 2016).
%
%   t = RUN.ML_TSCAN(X) returns the pseudotime of every cell, scaled to
%   [0, 1] along the main path of the minimum spanning tree over cluster
%   centres. X is genes-by-cells. Cells in clusters that do not lie on that
%   path get NaN: the method orders cells along one backbone, and a cell on
%   a side branch has no position on it.
%
%   T IS INDEXED BY CELL, matching PKG.I_PSEUDOTIME_BY_SPLINEFIT, the other
%   method SC_TRAJECTORY offers. It used to return the output of
%   [~, t] = sortrows(...), which is a sort index -- "which cell comes
%   k-th" -- not "cell k's pseudotime". CLI.CMD_TRAJECTORY writes that
%   column straight into a table beside the cell ids and stores it as the
%   'pseudotime' cell attribute, so every cell was being labelled with an
%   unrelated cell's index.

% ref: https://academic.oup.com/nar/article/44/13/e117/2457590
% load example_data\example10xdata.mat
% [X,genelist]=sc_selectg(X,genelist,7,8);
% [X]=sc_selectc(X);
% [t_pseudotime]=sc_tscan(X);

p = inputParser;
addRequired(p, 'X', @isnumeric);
addOptional(p, 'do_geneculst', true, @islogical);
addOptional(p, 'do_reduce', true, @islogical);
addOptional(p, 'plotit', true, @islogical);
parse(p, X, varargin{:});

do_geneculst = p.Results.do_geneculst;
do_reduce = p.Results.do_reduce;
plotit = p.Results.plotit;

if do_geneculst
    % in order to alleviate the effect of drop-out events (20) on
    % the subsequent analyses, genes with similar expression patterns are
    % grouped into clusters by hierarchical clustering (using Euclidean
    % distance and complete linkage). The number of clusters is set to
    % be 5% of the total number of genes with non-zero expression. For each
    % cluster and each cell, the expression measurements of all genes in
    % the cluster are averaged to produce a cluster-level expression which
    % will be used for subsequent MST construction.
    n = round(0.05*size(X, 1));
    Y = pdist(X);
    % 'complete' is what the paragraph above specifies and what TSCAN uses.
    % LINKAGE with no method is SINGLE linkage, which chains: it produces
    % one giant gene cluster plus a spray of singletons, so the averaging
    % below collapses most of the data into a single near-meaningless row
    % and the PCA that follows runs on the leftovers.
    Z = linkage(Y, 'complete');
    clu = cluster(Z, 'maxclust', n);
    % X = grpstats(X, clu, @(x) mean(x, 1));
    X = splitapply(@(x) mean(x, 1), X, clu);
    % Xn=[];
    % for k=1:27
    %    Xn=[Xn;mean(X(idx==k,:),1)];
    % end
end


if do_reduce
    % Briefly, Ei from all cells are organized into a H × N matrix E?. Each row
    % corresponds to a gene cluster. The matrix is standardized such that
    % expression values within each row have zero mean and unit standard deviation
    pcadim = 5;
    Xn = zscore(X, 0, 2);
    [~, X] = pca(Xn', 'NumComponents', pcadim);
end


% z A matrix whose [i,k]th entry is the probability that observation i in the test data
% belongs to the kth class. https://github.com/zji90/TSCAN/blob/master/R/exprmclust.R
% G The optimal number of mixture components
% There are other options you can use to help select the appropriate number of components for a Gaussian mixture model. For example,
% Compare multiple models with varying numbers of components using information criteria, e.g., AIC or BIC.
% Estimate the number of clusters using evalclusters, which supports, the Calinski-Harabasz criterion and the gap statistic, or other criteria.
% res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, modelNames = modelNames)
% https://www.mathworks.com/help/stats/clustering-using-gaussian-mixture-models.html

% RNG(271) is global state: leaving it set changes every later random draw
% in the session. ONCLEANUP puts the caller's stream back, including if
% anything below throws.
rngState = rng(271);
restoreRng = onCleanup(@() rng(rngState));

clunum = i_selectclunum(X);


% Use the same regularization as the AIC selection loop above; refitting
% without it can fail on an ill-conditioned covariance for the very
% clunum that selection just chose.
gmfit = fitgmdist(X, clunum, 'CovarianceType', 'full', 'RegularizationValue', 0.1);
clusterid = cluster(gmfit, X);
% clucenter = grpstats(X, clusterid, @mean);
clucenter = splitapply(@mean, X, clusterid);

txtc = strings(size(clucenter, 1), 1);
for k = 1:size(clucenter, 1)
    txtc(k) = string(sprintf('Clu%d', k));
end

G = graph(squareform(pdist(clucenter)), txtc);
T = minspantree(G);
D = distances(T);
[i, j] = find(D == max(D(:)));
clupath = shortestpath(T, i(1), j(1));

%%
% Order the cells of each cluster on the path along the direction the path
% travels through that cluster.
%
% The path used to be closed back to its own start, clupath(end+1) =
% clupath(1), so that the loop -- which reads clupath(k) and clupath(k+1)
% -- would give the terminal cluster a turn as k. It did, but with
% clupath(k+1) then being the ROOT, so the direction vector for that
% cluster pointed back down the trajectory and every cell in the final
% stretch came out in reverse. The terminal cluster now takes its direction
% from its predecessor, which is the same forward direction.
numSeg = numel(clupath);
tt = nan(size(X, 1), 2);
for k = 1:numSeg
    i = clupath(k);
    if k < numSeg
        difvec = clucenter(clupath(k+1), 1:3) - clucenter(i, 1:3);
    else
        difvec = clucenter(i, 1:3) - clucenter(clupath(k-1), 1:3);
    end

    idx = clusterid == i;
    if ~any(idx), continue; end

    difv = difvec / norm(difvec);
    projection = difv*X(idx, 1:3).';

    % The rank of each cell, not the sort index. SORT returns idxv with
    % idxv(r) = the cell of rank r; assigning that straight back labels
    % each cell with the identity of a different one -- the inverse
    % permutation. For projections [0.5 0.1 0.9 0.3] the old line stored
    % [2 4 1 3] where the ranks are [3 1 4 2].
    [~, ord] = sort(projection);
    withinRank = zeros(size(ord));
    withinRank(ord) = 1:numel(ord);

    tt(idx, 1) = k;
    tt(idx, 2) = withinRank;
end

% Cells in clusters off the main path have no position on it. They used to
% keep the NaN they were preallocated with and then be swept to the end by
% SORTROWS, which put them at the highest pseudotime of all rather than
% reporting them as unplaced -- so a user looking for late-trajectory genes
% got the side branch.
onPath = ~isnan(tt(:, 1));
if ~all(onPath)
    warning('ml_TSCAN:cellsOffMainPath', ...
        ['%d of %d cells lie in clusters off the main path and have no ', ...
        'pseudotime; they are returned as NaN.'], ...
        nnz(~onPath), numel(onPath));
end

t = nan(size(X, 1), 1);
placed = find(onPath);
[~, order] = sortrows(tt(placed, :), [1, 2]);
ordered = placed(order);
if numel(ordered) > 1
    t(ordered) = (0:numel(ordered)-1).'/(numel(ordered) - 1);
elseif isscalar(ordered)
    t(ordered) = 0;
end

if plotit
    subplot(2, 2, 2)
    p = plot(G); % ,'EdgeLabel',G.Edges.Weight);
    highlight(p, T, 'EdgeColor', 'r', 'LineWidth', 4.5)
    title('Graph of cell cluster centers')

    subplot(2, 2, 3)
    plot(T);
    title('MST')

    subplot(2, 2, 4)
    p = plot(T);
    highlight(p, clupath, 'EdgeColor', 'r', 'LineWidth', 4.5)
    title('Main path in MST')

    subplot(2, 2, 1)
    gui.i_gscatter3(X(:, 1:3), clusterid);
    hold on
    for k = 1:clunum
        plot3(clucenter(k, 1), clucenter(k, 2), clucenter(k, 3), '+', 'markersize', 20)
        text(clucenter(k, 1), clucenter(k, 2), clucenter(k, 3), sprintf('%d', k), 'fontsize', 30)
    end
    title('Cell clusters')
end

end


function clunum = i_selectclunum(X)
% Choose the number of mixture components by AIC over 2..9.
%
% Fitting eight mixtures produces a stream of convergence warnings that are
% expected and that the caller cannot act on, so they are silenced -- but
% only around this sweep. The original silenced them with a bare
% "warning off" in the main body and re-enabled with "warning on", which
% both left them off if anything in between threw, and would have swallowed
% the off-main-path warning further down. ONCLEANUP restores the caller's
% exact warning state on every path out.
warnState = warning('off', 'all');
restoreWarn = onCleanup(@() warning(warnState));

vaic = zeros(1, 8);
for k = 1:8
    try
        gm = fitgmdist(X, k+1, 'RegularizationValue', 0.1);
        vaic(k) = gm.AIC;
    catch
        vaic(k) = nan;
    end
end
[~, idx1] = min(vaic);
clunum = idx1 + 1;
end

%{
c1=clucenter(1,1:2);
c2=clucenter(2,1:2);
x1=X(clusterid==1,1:2);
difvec=c2-c1;
difv=difvec/norm(difvec);
[~,idx]=sort(difv*x1');
x1=x1(idx,:);
for k=1:size(x1,1)
    text(x1(k,1),x1(k,2),sprintf('%d',k));
end
%}
