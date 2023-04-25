function [es_max] = sctype(X, scaled, gs, gs2, gene_names_to_uppercase, varargin)


% marker sensitivity
[~, marker_stat] = sort(histcounts(gs), 'descend');
marker_sensitivity = table(scales.rescale(double(marker_stat), ...
    [0 1], [length(gs) 1]), ...
string(names(marker_stat)'), ...
'VariableNames', {'score_marker_sensitivity', 'gene_'});

% convert gene names to Uppercase
if gene_names_to_uppercase
    X.Properties.RowNames = upper(X.Properties.RowNames);
end

% subselect genes only found in data
names_gs_cp = string(names(gs)); 
names_gs_2_cp = string(names(gs2));
gs = cellfun(@(d_) {X.Properties.RowNames(ismember(string(X.Properties.RowNames), ...
    string(d_)))}, gs);
gs2 = cellfun(@(d_) {X.Properties.RowNames(ismember(string(X.Properties.RowNames), string(d_)))}, gs2);
gs = cell2struct(gs, names_gs_cp, 2); gs2 = cell2struct(gs2, names_gs_2_cp, 2);
cell_markers_genes_score = marker_sensitivity(ismember(marker_sensitivity.gene_, unique([gs{:} gs2{:}])), :);

% z-scale if not
if ~scaled
    Z = transpose(zscore(transpose(X))); 
else
    Z = X;
end

% multiple by marker sensitivity
for jj = 1:height(cell_markers_genes_score)
Z(string(cell_markers_genes_score(jj, 'gene_')), :) = ...
Z(string(cell_markers_genes_score(jj, 'gene_')), :) * cell_markers_genes_score(jj, 'score_marker_sensitivity');
end

% subselect only with marker genes
Z = Z(unique([gs{:} gs2{:}]), :);

% combine scores
es = cellfun(@(gss_) cellfun(@(j) ...
sum(Z(gss_, j), 'omitnan') / ...
sqrt(length(gss_)) + sum(Z(gs2.(gss_), j), ...
'omitnan') / sqrt(length(gs2.(gss_))), ...
num2cell(1:size(Z, 2))), names_gs_cp, ...
'UniformOutput', false);

es = vertcat(es{:});
es_max = array2table(es, 'RowNames', names_gs_cp, ...
    'VariableNames', X.Properties.VariableNames);
es_max = es_max(~all(ismissing(es_max), 2), :);

end