function MI_mat = MI_construction(data)
    % MI_construction computes the mutual information from data
    % INPUT:
    % data =====> Contains X count matrix and y target as follows
    %             data = [X; y]; (gene by (cell basis + target) )
    % USAGE:
    % R0 = MI_construction(data);
    % save('R0.mat', 'R0', '-v7.3');

    % Data is count matrix with genes in rows and cells in columns
    data = sparse(data);

    % transpose for efficient access (vectorization) across observations
    data = transpose(data);
    nobs = size(data,1);
    ngene = size(data,2); 
    MI_mat = zeros(ngene);
    % MI across data's rows (Computing upper triangular) in parallel 
    % NOTE: Diagonal elements contain zeros
    parfor jg = 1:ngene
        tmp = zeros(ngene,1);
        for ig = jg + 1: ngene
            % Computing pair MI with binning distribution
            tmp(ig) = BinPairMI( full(data(1:nobs,jg)), ...
                                 full(data(1:nobs,ig)) );
        end
        MI_mat(jg,:) = tmp;
    end
    % Copy upper triangular to lower triangular
    MI_mat = MI_mat + triu(MI_mat, 1)';
end