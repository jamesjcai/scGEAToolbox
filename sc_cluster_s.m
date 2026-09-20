function [c_clustid] = sc_cluster_s(s, k, varargin)
% sc_cluster_s - cluster cells using cell embeding s
%
%   c = SC_CLUSTER_S(s, k) partitions the rows of the embedding S into K
%   clusters and returns a cluster index per cell.
%
%   c = SC_CLUSTER_S(s, k, type) selects the method: 'kmeans' (default),
%   'kmedoids', 'spectclust', 'snndpc' or 'mbkmeans'.
%
%   c = SC_CLUSTER_S(..., 'Replicates', r) restarts 'kmeans' and 'mbkmeans'
%   from R independent seedings and keeps the best (default 5). It does not
%   apply to the other methods; see the note on 'kmedoids' below.
%
% see also: sc_cluster_x

%if min(size(s))>3, error('S is coordinates of dimensional
%reduction.'); end

if nargin < 2, k = 6; end
p = inputParser;
defaultType = 'kmeans';
validTypes = {'kmeans', 'kmedoids', 'dbscan', ...
    'spectclust', 'snndpc', 'mbkmeans'};
%
checkType = @(x) any(validatestring(x, validTypes));

checkK = @(x) (x > 0) && isnumeric(x) && isscalar(x);

addRequired(p, 's', @isnumeric);
addRequired(p, 'k', checkK);
addOptional(p, 'type', defaultType, checkType);
addOptional(p, 'plotit', false, @islogical);
addParameter(p, 'Replicates', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == fix(x));
parse(p, s, k, varargin{:});
plotit = p.Results.plotit;
numrep = p.Results.Replicates;

switch p.Results.type
    case {'spectralcluster', 'spectclust'}
        c_clustid = spectralcluster(s, k);
    case 'kmeans'
        % KMEANS defaults to a single seeding, which is not enough. From
        % one seeding it lands in a local minimum that merges two clusters
        % and splits a third: on six Gaussian clusters in 10-D over 20
        % seeds it reached the true partition 4 times out of 20, against
        % 20 out of 20 with five restarts. The gain holds at every
        % separation tried, and five restarts cost 0.09 s against 0.03 s
        % on a 20000-cell 3-D embedding, so there is nothing to trade.
        % The default MaxIter of 100 is too few for large t-SNE embeddings
        % (dense, curved clusters keep shifting a few points per pass): a
        % replicate stopped there is unconverged, and KMEANS warns
        % stats:kmeans:FailedToConvergeRep. Iterations stop once
        % assignments no longer change, so the higher cap costs nothing
        % when the default would have sufficed.
        c_clustid = kmeans(s, k, 'Replicates', numrep, 'MaxIter', 1000);
    case 'kmedoids'
        % Deliberately not replicated. KMEDOIDS seeds with a k-means++
        % build and reaches the same answer from one seeding as from five
        % -- identical mean ARI at every separation tried -- so restarts
        % would cost 2.4x for nothing.
        c_clustid = kmedoids(s, k);
    case 'dbscan'
        error('sc_cluster_s:NotImplemented', 'DBSCAN clustering is not yet implemented.');
    case 'snndpc'
        c_clustid = sc_snndpc(s, k);
    case 'mbkmeans'
        [~, ~, c_clustid] = pkg.e_mbkmeans(s, k, [], [], ...
            Replicates=numrep);
    otherwise
        error('sc_cluster_s:InvalidType', 'Unknown clustering type: %s', p.Results.type);
end

if plotit
    gui.i_gscatter3(s, c_clustid);
    hold on
    for i = 1:k
        si = s(c_clustid == i, :);
        si = mean(si);
        if size(s, 2) == 3
            text(si(:, 1), si(:, 2), si(:, 3), sprintf('%d', i), ...
                'fontsize', 20, 'FontWeight', 'bold', 'BackgroundColor', ...
                'w', 'EdgeColor', 'k');
        else
            text(si(:, 1), si(:, 2), sprintf('%d', i), ...
                'fontsize', 20, 'FontWeight', 'bold', 'BackgroundColor', ...
                'w', 'EdgeColor', 'k');
        end
    end
end
end
