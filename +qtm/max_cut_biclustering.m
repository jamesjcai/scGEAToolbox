function [gene_labels, sample_labels] = max_cut_biclustering(expression_matrix, k)
    % Initial clustering of genes
    [gene_labels, ~] = kmeans(expression_matrix, 2);
    gene_labels = gene_labels - 1;  % Convert to 0 and 1
    
    % Iterative refinement
    for i = 1:k
        % Update sample clustering based on gene clusters
        sample_scores = sum(expression_matrix .* (2*gene_labels - 1), 1);
        sample_labels = double(sample_scores > 0);
        
        % Update gene clustering based on sample clusters
        gene_scores = sum(expression_matrix .* (2*sample_labels - 1), 2);
        gene_labels = double(gene_scores > 0);
    end
end

%{ 
% Example usage
expression_matrix = rand(100, 20);  % 100 genes, 20 samples
[gene_clusters, sample_clusters] = qtm.max_cut_biclustering(expression_matrix, 10);

disp('Gene clusters:');
disp(gene_clusters');
disp('Sample clusters:');
disp(sample_clusters);
%}