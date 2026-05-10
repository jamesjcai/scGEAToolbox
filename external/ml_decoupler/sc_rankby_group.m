function T = sc_rankby_group(scores, terms, labels, options)
% SC_RANKBY_GROUP  Identify top-enriched gene sets per cell cluster.
%
% For each cell group (cluster / cell type), ranks gene sets by their
% mean enrichment score and returns the top N per group. Equivalent to
% decoupler-py dc.tl.rankby_group().
%
% USAGE:
%   T = sc_rankby_group(scores, terms, labels)
%   T = sc_rankby_group(scores, terms, labels, TopN=5)
%
% INPUTS:
%   scores - n_cells x n_terms enrichment score matrix (from sc_ulm)
%   terms  - 1 x n_terms string array of term names (from sc_ulm)
%   labels - n_cells x 1 vector of group labels (numeric, string, or categorical)
%
% OUTPUT:
%   T - table with columns:
%       Group       - group/cluster label
%       Rank        - rank within group (1 = best)
%       Term        - gene set / pathway / cell type name
%       MeanScore   - mean ULM t-statistic in this group
%       PctPositive - fraction of cells in group with score > 0
%       MeanVsRest  - mean score in group minus mean score in all other cells
%
% OPTIONS:
%   TopN (default 3) - number of top terms to report per group
%
% EXAMPLE:
%   net = sc_net_progeny("human");
%   [scores, ~, terms] = sc_ulm(sce, net);
%   T = sc_rankby_group(scores, terms, sce.c_cluster_id, TopN=5);
%   disp(T(T.Group == "0", :))

arguments
    scores  (:,:) double
    terms   (1,:) string
    labels  (:,1)
    options.TopN (1,1) double {mustBePositive, mustBeInteger} = 3
end

[n_cells, n_terms] = size(scores);
assert(numel(terms) == n_terms, ...
    'sc_rankby_group: length of terms (%d) must match columns of scores (%d).', ...
    numel(terms), n_terms);
assert(numel(labels) == n_cells, ...
    'sc_rankby_group: length of labels (%d) must match rows of scores (%d).', ...
    numel(labels), n_cells);

labels  = categorical(string(labels(:)));
groups  = categories(labels);
n_groups = numel(groups);

rows = cell(n_groups, 1);
for g = 1:n_groups
    inGroup  = labels == groups{g};
    gScores  = scores(inGroup, :);          % n_in x n_terms
    rScores  = scores(~inGroup, :);         % n_out x n_terms

    meanIn   = mean(gScores, 1, 'omitnan');
    meanRest = mean(rScores, 1, 'omitnan');
    pctPos   = mean(gScores > 0, 1, 'omitnan');
    vsRest   = meanIn - meanRest;

    % Rank by mean score descending
    [~, order] = sort(meanIn, 'descend');
    topN = min(options.TopN, n_terms);
    topIdx = order(1:topN);

    groupRows = table( ...
        repmat(string(groups{g}), topN, 1), ...
        (1:topN)', ...
        terms(topIdx)', ...
        meanIn(topIdx)', ...
        pctPos(topIdx)', ...
        vsRest(topIdx)', ...
        'VariableNames', {'Group','Rank','Term','MeanScore','PctPositive','MeanVsRest'});
    rows{g} = groupRows;
end

T = vertcat(rows{:});

end
