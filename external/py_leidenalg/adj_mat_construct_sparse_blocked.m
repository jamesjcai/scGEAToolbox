function adjX = adj_mat_construct_sparse_blocked(sce, method, K, chunk_size)
    % INPUT:
    % sce --------> SCE object 
    % method -----> Neighbor method ('knn' or 'mnn')
    % K ----------> Number of neighbors
    % chunk_size -> Number of rows to process in each chunk
    % OUTPUT:
    % adjX -------> Sparse adjacency matrix
    % AUTHOR: Selim Romero, Texas A&M University

    % Input validation
    method = lower(method);
    if ~ismember(method, {'knn', 'mnn'})
        error('Method must be either ''knn'' or ''mnn''.');
    end
    if nargin < 3; K = 15; end
    if ~isnumeric(K) || K <= 0 || K ~= round(K)
        error('K must be a positive integer.');
    end
    if nargin < 4; chunk_size = 1e4; end
    if ~isnumeric(chunk_size) || chunk_size <= 0 || chunk_size ~= round(chunk_size)
        error('chunk_size must be a positive integer.');
    end

    fprintf("Computing chunked adjacency matrix! \n");
    tic;
    % Normalize and preprocess input data (genes by cells)
    X = sce.X; 
    X = sc_norm(X, 'type', 'libsize'); 
    X = log1p(X)'; % (cells by genes)
    [U, ~, ~] = svds(X', 50); % (cells bu meta-genes)
    X = X * U;
    fprintf("Time for preparing data: %f \n", toc);

    % Initialize sparse adjacency matrix
    num_cells = size(X, 1);
    adjX = sparse(num_cells, num_cells);

    tic;
    % Chunk-wise computation of pairwise similarities
    num_chunks = ceil(num_cells / chunk_size);
    for chunk_idx = 1:num_chunks
        fprintf("Processing chunk %d of %d...\n", chunk_idx, num_chunks);

        % Define the range of rows for the current chunk
        ibeg = (chunk_idx - 1) * chunk_size + 1;
        iend = min(chunk_idx * chunk_size, num_cells);

        % Compute similarities for the current chunk
        % Cosine similarity
        %sim_chunk = 1 - pdist2(X(ibeg:iend, :), X, 'cosine');
        % Euclidean distance
        sim_chunk = 1 - pdist2(X(ibeg:iend, :), X, 'euclidean');

        switch method
            case 'knn'
                % KNN for the current chunk
                for i = 1:size(sim_chunk, 1)
                    global_idx = ibeg + i - 1;
                    [~, idx] = sort(sim_chunk(i, :), 'descend');
                    adjX(global_idx, idx(1:K)) = 1;
                end
            case 'mnn'
                % MNN for the current chunk
                for i = 1:size(sim_chunk, 1)
                    global_idx = ibeg + i - 1;
                    [~, idx] = sort(sim_chunk(i, :), 'descend');
                    neighbors_i = idx(1:K);
                    for j = neighbors_i
                        [~, idx_j] = sort(sim_chunk(:, j), 'descend');
                        neighbors_j = idx_j(1:K);
                        if ismember(global_idx, neighbors_j)
                            adjX(global_idx, j) = 1;
                            adjX(j, global_idx) = 1;
                        end
                    end
                end
        end
    end

    % Ensure the adjacency matrix is symmetric
    adjX = max(adjX, adjX');
    fprintf("Total time for neighbors: %f \n", toc);

    % Final adjacency matrix is sparse
    adjX = sparse(adjX);
end
