function [X, genelist] = sc_rmdugenes(X, genelist, methodid)

if nargin<3, methodid = 1; end

    % Remove genes with duplicate name - Optimized version
    [genelist_out, first_idx, group_idx] = unique(genelist, 'stable');

    num_duplicates = length(genelist) - length(first_idx);
    if num_duplicates == 0
        return;
    else
        switch methodid
            case 1
                genelist_out = matlab.lang.makeUniqueStrings(genelist);
            case 2
                X = X(first_idx, :);                
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
            case 4
                tic;
                X = splitapply(@(rows) sum(rows,1), X, group_idx);
                toc;
            case 5
                tic;
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
                toc;
        end
        genelist = genelist_out;        % genelist(first_idx);        
        % warning('Duplicate gene names are removed. %d duplicate genes found.', num_duplicates);
    end

end

%{
 Pre-allocate result matrix
X_new = zeros(length(first_idx), size(X, 2));

% Vectorized summation for each unique gene group
for i = 1:length(first_idx)
    group_members = (first_idx == i);
    X_new(i, :) = sum(X(group_members, :), 1);
end

% Update outputs
X = X_new;
genelist = genelist(first_idx);
%}



%{
function [X, genelist] = sc_rmdugenes(X, genelist)
% Remove genes with duplicate name

[~, w] = unique(genelist, 'stable');
duplicate_indices = setdiff(1:numel(genelist), w);
if isempty(duplicate_indices)
    return;
end

for k = 1:length(duplicate_indices)
    idx = find(genelist == genelist(duplicate_indices(k)));
    X(idx(1), :) = sum(X(idx, :), 1);
end
X(duplicate_indices, :) = [];
genelist(duplicate_indices) = [];
end
%}

%{
function [X, genelist] = sc_rmdugenes(X, genelist)
% Remove genes with duplicate name - Optimized version
[~, first_idx, group_idx] = unique(genelist, 'stable');

% If no duplicates, return early
if length(first_idx) == length(genelist)
    return;
end

% Pre-allocate result matrix
X_new = zeros(length(first_idx), size(X, 2));

% Vectorized summation for each unique gene group
for i = 1:length(first_idx)
    group_members = (group_idx == i);
    if sum(group_members) == 1
        X_new(i, :) = X(group_members, :);
    else
        X_new(i, :) = sum(X(group_members, :), 1);
    end
end

% Update outputs
X = X_new;
genelist = genelist(first_idx);
end
%}

%{

function [X_out, genelist_out] = sc_rmdugenes(X, genelist)
% SC_RMDUGENES  Remove duplicate gene names and sum their expression rows.
% Uses splitapply when available (R2020b+), otherwise falls back to loop.
%
% Syntax:
%   [X_out, genelist_out] = sc_rmdugenes(X, genelist)
%
% Inputs:
%   X         - Numeric gene expression matrix (genes × cells)
%   genelist  - Cell array of gene names (genes × 1)
%
% Outputs:
%   X_out        - Matrix with duplicates summed
%   genelist_out - Unique gene list in order of first appearance

    arguments
        X {mustBeNumeric, mustBeNonempty}
        genelist {mustBeVector}
    end

    % Map each gene to its group index
    [genelist_out, ~, ic] = unique(genelist, 'stable');

    % Check if splitapply exists and supports matrix data
    if exist('splitapply', 'file') == 2
        try
            % Modern vectorized approach
            X_out = splitapply(@(rows) sum(rows, 1), X, ic);
            return
        catch
            % If splitapply fails for older versions or data types, fall back
        end
    end

    % Fallback loop (fast because it iterates over unique genes)
    X_out = zeros(numel(genelist_out), size(X, 2), class(X));
    for k = 1:numel(genelist_out)
        X_out(k, :) = sum(X(ic == k, :), 1);
    end
end

%}