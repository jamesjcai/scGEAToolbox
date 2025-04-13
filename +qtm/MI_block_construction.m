function MI_mat = MI_block_construction(data, data2, comp_ii )
    % MI_construction computes the mutual information from data
    % INPUT:
    % data =======> Contains X count matrix and y target as follows
    %               data = [X; y]; (gene by (cell basis + target) )
    % data2 ======> (Optional)
    %               Contains X2 count matrix and y target as follows
    %               data = [X2; y2]; (gene by cell basis + target by cell basis) )
    % comp_ii ====> Compute the ig-ig diagonal terms? true/false
    % OUTPUT:
    % MI_mat =====> Mutual 
    % USAGE:
    % R0 = MI_construction(data);
    % save('R0.mat', 'R0', '-v7.3');
    
    % Data is count matrix with genes in rows and cells in columns
    data = sparse(data);

    if nargin < 2 || isempty(data2)
        data2 = data;
    else
        data2 = sparse(data2);
    end
    if nargin < 3 ||  isempty(comp_ii); comp_ii = false; end

    % transpose for efficient access (vectorization) across observations
    data = data';
    data2 = data2';

    nobs = size(data,1);
    nobs2 = size(data2,1);
    if nobs ~= nobs2
        error("Number of observations in MI_block_construnction do not match")
    end
    ngene = size(data, 2); 
    ngene2 = size(data2, 2); 

    nbatch = 500;
    nblock = ceil(ngene/ nbatch);
    nblock2 = ceil(ngene2/ nbatch);

    tic;
    MI_mat = zeros(ngene, ngene2);
    % MI across data's rows (Computing upper triangular) in parallel 
    % NOTE: Diagonal elements contain zeros
    progressbar('MI_block_construction begins: ');
    progressbar(0);
    for iblock = 1:nblock
        ibeg = 1 + (iblock - 1)*nbatch;
        iend = min( iblock*nbatch, ngene);
        % Cells + target by gene block
        data_block = data(:, ibeg:iend);
        %fprintf("Iblock %d , ibeg %d  and iend %d \n ", iblock, ibeg, iend);
        for jblock = iblock:nblock2

            jbeg = 1 + (jblock - 1)*nbatch;
            jend = min( jblock*nbatch, ngene2);
            %fprintf("Jblock %d , jbeg %d  and jend %d \n ", jblock, jbeg, jend);

            % Cells + target by gene block2
            data_block2 = data2(:, jbeg:jend);
            if iblock == jblock
                mode = 1;
            else
                mode = 2; 
            end
            MI_mat(ibeg:iend, jbeg:jend) = ...
                MI_block( data_block, data_block2, mode, comp_ii );
    
        end
        progressbar(iblock/nblock*100);
    end
    mi_time = toc;
    fprintf("MI time %f \n", mi_time);

    % Copy upper triangular to lower triangular if it is square matrix
    if ~comp_ii && ngene2 == ngene
        MI_mat = MI_mat + triu(MI_mat, 1)';
    end
end