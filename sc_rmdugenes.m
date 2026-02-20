function [X, genelist] = sc_rmdugenes(X, genelist, methodid)

if nargin<3, methodid = 1; end

    % Remove genes with duplicate name - Optimized version
    [genelist_out, first_idx, group_idx] = unique(genelist, 'stable');

    num_duplicates = length(genelist) - length(first_idx);
    if num_duplicates == 0
        return;
    else
        fprintf('Duplicate genes (n = %d) found.\n', num_duplicates);
        switch methodid
            case 1
                genelist_out = matlab.lang.makeUniqueStrings(genelist);
                disp('Duplicate genes are renamed by appending an underscore and a number to duplicates.')
            case 2
                X = X(first_idx, :);
                disp('Duplicate genes are removed; only the first occurrence of duplicates are kept.')
            case 3

                nGenes = numel(first_idx);
                nCells = size(X, 2);
                
                if issparse(X)
                    % Extract nonzero elements
                    [row, col, val] = find(X);
                    % Map each original row index to its unique-gene index
                    row = group_idx(row);
                    % Sum duplicates directly into sparse matrix
                    if isMATLABReleaseOlderThan('R2025a')
                        X = sparse(row, col, val, nGenes, nCells);
                    else
                        X = single(sparse(row, col, val, nGenes, nCells)); % lower memory in R2025a+
                    end                    
                else
                    % Dense case: loop over columns, group sum
                    X_new = zeros(nGenes, nCells, class(X));
                    for j = 1:nCells
                        X_new(:, j) = accumarray(group_idx, X(:, j), [nGenes, 1]);
                    end
                    X = X_new;
                end                 
                % GPT5 solution! 100x faster
                disp('Duplicate genes are collapsed/merged by summing their expression rows')
            case 4                
                X = splitapply(@(rows) sum(rows,1), X, group_idx);                
            case 5                
                X_new = zeros(length(first_idx), size(X, 2), class(X));
                % consider grplen = groupcounts(group_idx);
                for ix = 1:length(first_idx)
                    group_members = (group_idx == ix);
                    if sum(group_members) == 1
                        X_new(ix, :) = X(group_members, :);
                    else
                        X_new(ix, :) = sum(X(group_members, :), 1);
                    end
                end
                if issparse(X), X_new = sparse(X_new); end
                X = X_new;                
        end
        genelist = genelist_out;        % genelist(first_idx);        
        % warning('Duplicate gene names are removed. %d duplicate genes found.', num_duplicates);
    end

end
