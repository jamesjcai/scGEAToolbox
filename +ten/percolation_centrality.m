function perc_centrality = percolation_centrality(adj_matrix, percolation_prob)
    % adj_matrix: adjacency matrix of the graph
    % percolation_prob: vector of percolation probabilities for each node

    % Number of nodes
    n = size(adj_matrix, 1);
    
    % Initialize percolation centrality vector
    perc_centrality = zeros(n, 1);
    
    % Iterate over each node to calculate its percolation centrality
    for i = 1:n
        % Initialize a queue for BFS
        queue = i;
        visited = false(n, 1);
        visited(i) = true;
        
        % BFS to calculate the reachability considering percolation probabilities
        while ~isempty(queue)
            current_node = queue(1);
            queue(1) = [];
            
            neighbors = find(adj_matrix(current_node, :));
            for neighbor = neighbors
                if ~visited(neighbor) && rand() < percolation_prob(neighbor)
                    visited(neighbor) = true;
                    queue(end + 1) = neighbor; %#ok<AGROW>
                end
            end
        end
        
        % Percolation centrality is the sum of reachable nodes
        perc_centrality(i) = sum(visited);
    end
end


%{
% Example adjacency matrix for a simple graph
adj_matrix = [0 1 0 0;
              1 0 1 1;
              0 1 0 1;
              0 1 1 0];

% Example percolation probabilities for each node
percolation_prob = [0.8, 0.6, 0.7, 0.9];

% Calculate percolation centrality
perc_centrality = ten.percolation_centrality(adj_matrix, percolation_prob);

% Display the results
disp('Percolation Centrality:');
disp(perc_centrality);
%}