function MI_mat = MI_construction(data)
    % MI_construction computes the mutual information from data
    % INPUT:
    % data =====> Contains X count matrix and y target as follows
    %             data = [X; y]; (gene by (cell basis + target) )
    % USAGE:
    % R0 = MI_construction(data);
    % save('R0.mat', 'R0', '-v7.3');

    % Data is count matrix with genes in rows and cells in columns
    % data = sparse(data);

    % transpose for efficient access (vectorization) across observations
    data = transpose(data);
    [ngenep1] = size(data, 2);    

    MI_mat = zeros(ngenep1);
    % MI across data's rows (Computing upper triangular) in parallel 
    % NOTE: Diagonal elements contain zeros
    for jg = 1:ngenep1
        tmp = zeros(ngenep1, 1);
        d_jg = data(:, jg);
        parfor ig = jg + 1: ngenep1
            tmp(ig) = qtm.BinPairMI( d_jg, data(:, ig) );
        end
        MI_mat(jg,:) = tmp;
    end
    MI_mat = MI_mat + triu(MI_mat, 1)';
end