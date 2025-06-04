function [lgu, dropr, lgcv, genes, X_sorted, ...
    removedgeneidx, removedT] = sc_genestat_sparse_blocked(X, genelist, ...
    sortit, removeinf, blockSize)
% SC_GENESTAT_SPARSE_BLOCKED Computes gene statistics for scRNA-seq data,
% optimized for large, sparse matrices by processing cells in blocks.
%
% INPUTS:
%   X          - Sparse matrix of gene expression (genes x cells).
%   genelist   - Cell array or string array of gene names.
%   sortit     - Logical, whether to sort genes based on lgu, lgcv, dropr.
%                (Default: true)
%   removeinf  - Logical, whether to remove genes with Inf/NaN values in
%                lgu or lgcv. (Default: true)
%   blockSize  - Integer, number of cells (columns) to process in each block.
%                (Default: 10000)
%
% OUTPUTS:
%   lgu        - log1p(mean expression) for each gene.
%   dropr      - Dropout rate for each gene.
%   lgcv       - log1p(coefficient of variation) for each gene.
%   genes      - Sorted or filtered gene list.
%   X_sorted   - Sorted or filtered X matrix (if sortit is true).
%   removedgeneidx - Indices of genes removed due to Inf/NaN values.
%   removedT   - Table containing information about removed genes.
if nargin < 5 || isempty(blockSize), blockSize = 10000; end
if nargin < 4, removeinf = true; end
if nargin < 3, sortit = true; end
[numGenes, numCells] = size(X);
if nargin < 2 || isempty(genelist)
    genelist = "Gene"+string(1:numGenes).';
end
genelist = genelist(:);
if length(genelist) ~= numGenes
    error('Length of genelist must match the number of rows in X.');
end
geneidx_orig = (1:numGenes).';
% --- Calculate dropout rate (dropr) in blocks ---
nnz_per_gene = zeros(numGenes, 1, 'double'); % Ensure double for potential large sums
for k = 1:blockSize:numCells
    blockEnd = min(k + blockSize - 1, numCells);
    nnz_per_gene = nnz_per_gene + sum(X(:, k:blockEnd) > 0, 2);
end
dropr = 1 - nnz_per_gene ./ numCells;
% --- Calculate mean (u) and standard deviation (std_val) in blocks ---
% For mean: sum_x_per_gene / numCells
% For std: sqrt( (sum_x_sq_per_gene - (sum_x_per_gene.^2)/numCells ) / (numCells -1) )
% To avoid creating large intermediate matrices for X.^2, we sum X and X.^2
% iteratively.
sum_x_per_gene = zeros(numGenes, 1, 'double');
sum_x_sq_per_gene = zeros(numGenes, 1, 'double');
for k = 1:blockSize:numCells
    blockEnd = min(k + blockSize - 1, numCells);
    X_block = X(:, k:blockEnd); % This is a sparse sub-matrix
    sum_x_per_gene = sum_x_per_gene + full(sum(X_block, 2)); % Convert sum to full for addition
    % For sum of squares, operate on non-zero elements to maintain sparsity benefits
    % as much as possible before summing.
    % If X_block is very wide, X_block.^2 could become dense if many elements are non-zero.
    % However, for scRNA-seq, it's usually very sparse.
    sum_x_sq_per_gene = sum_x_sq_per_gene + full(sum(X_block.^2, 2));
end
u = sum_x_per_gene ./ numCells;
% Calculate variance and then standard deviation
% Handle the case where numCells = 1 to avoid division by zero in variance
if numCells > 1
    variance = (sum_x_sq_per_gene - (sum_x_per_gene.^2)./numCells) ./ (numCells - 1);
else
    variance = zeros(numGenes, 1); % Or NaN, depending on desired behavior for single cell
end
% Ensure non-negativity for variance due to potential floating point inaccuracies
variance(variance < 0) = 0;
std_val = sqrt(variance);
cv = std_val ./ u;
cv(u == 0) = NaN; % Avoid division by zero if mean is zero; assign NaN
lgu = log1p(u);
lgcv = log1p(cv);
genes = genelist;
X_sorted = X; % Initialize X_sorted
removedgeneidx = []; % Initialize to empty double array
removedT = [];     % Initialize to empty table
if sortit
    % Create a temporary array for sorting, handling NaNs by moving them to the end
    % For sorting, typically we want NaNs last.
    % A common way is to replace NaNs with Inf for sorting purposes if sorting ascending
    % or with -Inf if sorting descending, or use a custom sort.
    % Here, we'll sort and then handle NaNs in the filtering step if removeinf is true.
    
    sortData = [lgu, dropr, lgcv];
    
    % Sort rows: 1st col (lgu) asc, 3rd col (lgcv) asc, 2nd col (dropr) asc
    % Note: sortrows sorts NaNs last by default for ascending order.
    [~, si] = sortrows(sortData, [1, 3, 2]);
    
    lgu = lgu(si);
    dropr = dropr(si);
    lgcv = lgcv(si);
    genes = genelist(si);
    X_sorted = X(si, :); % Sort the original sparse matrix X
    geneidx_orig = geneidx_orig(si); % Keep track of original indices after sorting
else
    X_sorted = X; % If not sorting, X_sorted is the original X
end
% --- Remove Inf/NaN values ---
% si_nan is relative to the current state of 'genes', 'lgu', etc. (possibly sorted)
si_nan = isnan(lgu) | isinf(lgu) | isnan(lgcv) | isinf(lgcv);
if removeinf && any(si_nan)
    % Capture details of removed genes
    % Note: The zeros and ones in removedT from the original function might
    % represent specific flags. Replicating that structure.
    % If these columns had specific meanings, they should be documented.
    removedT = table(genes(si_nan), lgu(si_nan), lgcv(si_nan), dropr(si_nan), ...
        zeros(sum(si_nan),1,'double'), ... % Assuming double, adjust if needed
        ones(sum(si_nan),1,'double'), ...
        ones(sum(si_nan),1,'double'),...
        zeros(sum(si_nan),1,'double'), ...
        'VariableNames', {'Gene', 'LGU_removed', 'LGCV_removed', 'DropoutRate_removed', ...
                          'FlagCol1', 'FlagCol2', 'FlagCol3', 'FlagCol4'}); % Provide names
    removedgeneidx = geneidx_orig(si_nan); % Original indices of the removed genes
    lgu(si_nan) = [];
    lgcv(si_nan) = [];
    dropr(si_nan) = [];
    genes(si_nan) = [];
    if ~isempty(X_sorted) % Check if X_sorted is empty before trying to index
       X_sorted(si_nan, :) = [];
    end
    % geneidx_orig is already filtered implicitly by not including si_nan here for output indices
else
    removedgeneidx = zeros(0,1,'double'); % Ensure it's a 0x1 double if nothing removed
    removedT = table('Size',[0 8], ...
        'VariableTypes', {'string','double','double','double','double','double','double','double'}, ...
        'VariableNames', {'Gene', 'LGU_removed', 'LGCV_removed', 'DropoutRate_removed', ...
                          'FlagCol1', 'FlagCol2', 'FlagCol3', 'FlagCol4'}); % Empty table with schema
end
% If not sorting, X was assigned to X_sorted. If also not removing inf,
% X_sorted remains the original X. If removing inf but not sorting,
% X_sorted will have rows removed.
% The output X is now X_sorted to reflect changes.
end
 