function [Tsol, xsol] = qfeatures_qubo_base(X, g, y, K, readr)
    % qfeatures_qubo_base computes the QUBO feature selection (FS)
    % from count matrix X, genes g and y target.
    % INPUT:
    % X =====> Single cell coung matrix
    % g =====> Genes/features from single cell experiment
    % y =====> Predictor/target variable
    % K =====> Number of features to retrieve
    % readr => Read mutual information (MI) matrix? true/false
    % OUTPUT: 
    % Tsol ==> MATLAB table containing features and computation time
    
    if nargin < 5, readr = false; end
    tic;
    if readr
        % Loading MI if pointed out
        load('R0.mat')
    else
        % Computing MI
        data = sparse( [X; y] );
        R0 = MI_construction(data);
        clear data X y;
        save('R0.mat', 'R0', '-v7.3');
    end
    time_mi = toc;
    fprintf("MI construction time: %f \n", time_mi);

    fprintf("Looking for %d genes \n", K);
    % Redundancy matrix
    R = R0(1:end-1, 1:end-1) / (K - 1);
    % Importance vector
    J = R0(end, 1:end-1);

    tic;
    % Find optimal alpha to balance R and J with QUBO solver
    fun = @(alpha)howmany(alpha, R, J) - K;
    try
        alphasol = fzero(fun, [0 1]);
        fprintf("Alpha value fzero: %f \n", alphasol);
        [~, xsol] = howmany(alphasol, R, J);
    catch
        alphasol = fminsearch(fun, 1);
        fprintf("Alpha value fminsearch: %f \n", alphasol);
        [~, xsol] = howmany(alphasol, R, J);
    end
    time_zerof = toc;

    % Obtain QUBO function value
    fval = xsol.BestFunctionValue;
    % Re-compute Q  (balance of R and J)
    Q = (1 - alphasol) * R - alphasol * diag(J);
    % Function value per feature accoring best solution and Q matrix
    global_rank = Q * xsol.BestX;
    % Obtain top most important features
    [~, sort_idx] = sort(global_rank, 'ascend');
    % Re-sort original features
    gtmp = g(sort_idx);
    jdx = xsol.BestX(sort_idx) == 1;
    sol_genes = gtmp(jdx);

    % Make a row for features (Tsol requires it)
    nx = size(sol_genes, 1);
    if nx > 1
        sol_genes = gtmp(jdx)';
    end

    fprintf("QUBO FS alpha solver time: %f \n", time_zerof);

    Tsol = table(sol_genes, time_mi, time_zerof, fval, alphasol);
end