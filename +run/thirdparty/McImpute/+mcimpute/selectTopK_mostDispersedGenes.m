function [featureSelected_data, topK_indices] = selectTopK_mostDispersedGenes(gene_filtered_normed_data, K)

X = gene_filtered_normed_data;
%find cv2 for all genes
cv2 = var(X)./mean(X);

%divide into bins based on quantiles of mean cv2
[bin_numbers, mean_bin] = discretize(mean(X), [-inf, quantile(mean(X), 0.1:0.05:1), inf]  , 'IncludedEdge','right');
%[bin_numbers, mean_bin] = discretize(mean(X), [-inf, quantile(mean(X), 0.2:0.2:1), inf]  , 'IncludedEdge','right');


%divide cv2 of genes into groups based on bins and for each group find
%median and median adbsolute dispersion
mystats = @(x)[median(x); 1.4826*mad(x,1)];
var_by_bin = splitapply( mystats ,cv2,bin_numbers)'; %1st col has median, 2nd has MAD, nrows=#bins

%represent each gene by the median and MAD of its bin, and find normalized
%dispersion
bin_disp_median = var_by_bin(bin_numbers,1);
bin_disp_mad = var_by_bin(bin_numbers,2);
normalized_dispersion = abs(cv2'-bin_disp_median)./bin_disp_mad;

%pick top k indices having max dispersion
[sortedValues,sortIndex] = sort(-normalized_dispersion);
topK_indices = sortIndex(1:K);

%topKgenes = normalized_dispersion(indices);
featureSelected_data = gene_filtered_normed_data(:,topK_indices);
end

