function gene_order = tsp_gene_ordering(expression_matrix)
    % Calculate pairwise distances between gene expression profiles
    distances = pdist(expression_matrix, 'euclidean');
    distance_matrix = squareform(distances);
    
    % Solve TSP using MATLAB's optim.tsp function
    % Note: This requires the Optimization Toolbox
    gene_order = optim.tsp('CustomCost', @(x,y) distance_matrix(x,y));
end

% https://www.mathworks.com/help/optim/ug/traveling-salesman-problem-based.html

%{
Example usage
rng(42);  % Set random seed for reproducibility
expression_matrix = rand(20, 100);  % 20 genes, 100 samples
gene_order = qtm.tsp_gene_ordering(expression_matrix);

disp('Optimal gene order:');
disp(gene_order);

% Reorder the expression matrix
reordered_matrix = expression_matrix(gene_order, :);

disp('Shape of reordered matrix:');
disp(size(reordered_matrix));
%}