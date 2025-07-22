function sce = leiden_clustering_ann(sce, res, species, method)
    % leiden_clustering: Computes Leiden clustering interfaced from Python and annotates if desired
    % with Mutual Nearest Neighbors (MNN) or K-Nearest Neighbors (KNN).
    % INPUT:
    % sce -------> SCE object
    % method ----> Method to find neighbors ('mnn' or 'knn')
    % OUTPUT:
    % sce ------> SCE object containing Leiden clusters and corresponding annotation
    % Usage:
    % sce = leiden_annotation_sparse(sce, 'knn', 'mouse')
    % If no annotation is wanted, use
    % sce = leiden_annotation_sparse(sce, 'knn', [])
    % AUTHOR: Selim Romero, Texas A&M University
    
    if nargin < 2; res = 2.0; end
    if nargin < 3; species = []; end
    if nargin < 4; method = 'knn'; end

    species = lower(species);
    method = lower(method);

    fprintf("WARNING: sce object should be prior QC if desired...\n");
    fprintf("Leiden annotation with method: %s\n", method);

    % Set the Python environment (Python 3.11)
    env_bin = 'C:\Local_install\miniconda3\envs\scanpy_env_311\python.exe';
    %env_bin = 'F:\Anaconda\envs\scanpy_env_311\python.exe';
    if ispc
        env_bin = strrep(env_bin, "\", "\\");
    end
    % Linux format
    env_bin = "/home/ssromerogon/anaconda3/envs/scanpy_env_311/bin/python3";
    
    % Initialize the Python environment
    pe = pyenv('Version', env_bin);
    if pe.Status ~= "Loaded"
        fprintf("Reinitializing Python environment...\n");
        pe = pyenv('Version', env_bin);
        py.exec('import sys');
    end
    disp(pyenv);

    % Decide strategy for clustering based on memory
    [~, availableRAM, ~] = check_memory();
    mem_by_float = 8 / 1e9; % Memory per double-precision float (in GB)
    estimate_mem = sce.NumCells^2 * mem_by_float; % Estimate memory required for adjacency matrix

    % Determine if chunked strategy is needed
    if estimate_mem >= 0.7 * availableRAM
        fprintf("Chunked strategy... \n")
        chunked_strategy = true;
    else
        fprintf("Full strategy... \n")
        chunked_strategy = false;
    end

    % Select clustering method
    switch method
        case 'mnn'
            n_neighbors = 100; % Larger neighbors -> fewer clusters
            if chunked_strategy
                adjX = adj_mat_construct_sparse_blocked(sce, 'mnn', n_neighbors, 10000);
            else
                adjX = adj_mat_construct_sparse(sce, 'mnn', n_neighbors);
            end

        case 'knn'
            n_neighbors = 15; % Smaller neighbors -> more clusters
            if chunked_strategy
                adjX = adj_mat_construct_sparse_blocked(sce, 'knn', n_neighbors, 10000);
            else
                adjX = adj_mat_construct_sparse(sce, 'knn', n_neighbors);
            end

        otherwise
            error('Unknown method: %s. Method should be either ''mnn'' or ''knn''.', method);
    end

    disp(size(adjX));

    % Save adjacency matrix to file
    adj_file = 'adjX.txt';
    [i, j, val] = find(adjX);
    writematrix([i, j, val], adj_file, 'Delimiter', 'tab');
    clear adjX i j val;

    % Path to Python script
    python_executable = env_bin;
    leiden_wd = which('leiden_clustering_ann');
    leiden_wd = erase(leiden_wd, 'leiden_clustering_ann.m');
    if ispc
        leiden_wd = strrep(leiden_wd, "\", "\\");
    end
    python_script = strcat(leiden_wd, 'run_leiden_sparse.py');

    % Call Python script
    res_str = num2str(res);
    system_command = sprintf('%s %s %s %s', python_executable, python_script, adj_file, res_str);

    [status, cmdout] = system(system_command);

    % Clean up
    if exist(adj_file, 'file')
        delete(adj_file);
    end    

    % Handle Python script errors
    if status ~= 0
        disp('Error running the Python script:');
        disp(cmdout);
        return;
    else
        clusters = jsondecode(fileread('clusters.json'));
        delete('clusters.json');
        disp("Parsing Leiden:");
        disp(cmdout);
        nclus = length(unique(clusters));
        fprintf('Number Leiden clusters: %d\n', nclus);
    end

    % Assign clustering results to the SCE object
    sce.c_cluster_id = clusters + 1;

    % Optional: Annotate cell types
    tic;
    if ~isempty(species)
        fprintf("Annotating species: %s\n\n", species);
        sce = sce.assigncelltype(species, false);
    end
    time_assign = toc;
    fprintf("Time for cell annotation and embedding: %f\n", time_assign);
end
