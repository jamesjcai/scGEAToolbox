function [b] = qubofs(X, y, k, R0)

    if nargin <4, R0 = []; end
    if nargin < 3, k = 20; end
    if ~isrow(y), y = y.'; end

    data = [X; y];
    if isempty(R0)
        R0 = qtm.MI_construction(data);       
        % R0 = MI_block_construction(data);
    end

    % fprintf("Looking for %d genes \n", k);
    % Redundancy matrix
    R = R0(1:end-1, 1:end-1) / (k - 1);
    % Importance vector
    J = R0(end, 1:end-1);

    % Find optimal alpha to balance R and J with QUBO solver
    fun = @(alpha)qtm.howmany(alpha, R, J) - k;
    try
        alphasol = fzero(fun, [0 1]);
        fprintf("Alpha value fzero: %f \n", alphasol);
        [~, xsol] = qtm.howmany(alphasol, R, J);
    catch
        alphasol = fminsearch(fun, 1);
        fprintf("Alpha value fminsearch: %f \n", alphasol);
        [~, xsol] = qtm.howmany(alphasol, R, J);
    end

    % Obtain QUBO function value
    %fval = xsol.BestFunctionValue;
    % Re-compute Q  (balance of R and J)
    % Q = (1 - alphasol) * R - alphasol * diag(J);
    % Function value per feature accoring best solution and Q matrix
    b = logical(xsol.BestX);
    % global_rank = Q * b;



    %{    
    mdl = fitrtree(X',Y,CrossVal="on",Holdout=0.2);
    kfoldLoss(mdl)    
    
    mdl2 = fitrtree(X(b,:)',y,CrossVal="on",Holdout=0.2);
    kfoldLoss(mdl2)
    %}

    %{
    Obtain top most important features
    [global_rank, sort_idx] = sort(global_rank, 'ascend');
    % Re-sort original features
    gtmp = g(sort_idx);
    featureIndices = xsol.BestX(sort_idx) == 1;
    % featureIndices = find( xsol.BestX(sort_idx) == 1;)
    selectedGenes = gtmp(featureIndices);
    featureIndices = ismember( g, selectedGenes);

    % Make a row for features (Tsol requires it)
    if size(selectedGenes, 1) > 1
        selectedGenes = selectedGenes';
    end
    if size(featureIndices, 1) > 1
        featureIndices = featureIndices';
    end
    if size(global_rank, 1) > 1
        global_rank = global_rank';
    end
    fprintf("QUBO FS alpha solver time: %f \n", time_zerof);

    T = table(selectedGenes, featureIndices, global_rank, time_mi, time_zerof, fval, alphasol);
    %}
end