function boxplot_marker(data, allgenes, marker, cluster_labs, No_cluster)
% Box plot for each gene along all clusters
% Input:
%   -- data: Single cell data matrix where rows are associated with genes
%   and columns associated with cells
%   -- allgenes: Gene annotation
%   -- marker: Cell array of marker genes used for plot
%   -- cluster_labs: Cell clustering labels from SoptSC
%   -- No_cluster: Number of clusters
%   -- folder: folder name where the results saved
%

cmap1 = jet;
mymap1 = cmap1(1:end, :);
ncolor = size(mymap1, 1);
mycolor = mymap1(1:round(ncolor./No_cluster):ncolor, :);

group = cell(size(cluster_labs));
for i = 1:length(cluster_labs)
    group{i} = ['C', num2str(cluster_labs(i))];
end


ia = zeros(1, length(marker));
for i = 1:length(marker)
    for j = 1:length(allgenes)
        if strcmpi(marker{i}, allgenes{j})
            ia(i) = j;
        end
    end
end


gname = allgenes(ia);

display(allgenes(ia));

MM = data(ia, :);

for i = 1:No_cluster
    cluster_notation{i} = ['C', num2str(i)];
end

n = length(ia);
for i = 1:n
    figure;
    boxplot(MM(i, :), group, 'GroupOrder', cluster_notation, 'Notch', 'on', 'PlotStyle', 'compact', 'Widths', 0.9, 'Colors', mycolor, ...
        'LabelOrientation', 'horizontal');
    title(gname{i});
    % print([folder '\boxplot_mk_' gname{i}],'-dpdf','-r300');
end