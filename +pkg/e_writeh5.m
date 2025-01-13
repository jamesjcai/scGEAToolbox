function e_writeh5(X, genelist, filename, celltype, batchid)

% This function writes a sparse or dense matrix (X) along with optional metadata 
% (genelist, celltype, and s) into an HDF5 file.
% Parameters:
%   X         - Input matrix (can be sparse or dense)
%   genelist  - List of gene names (string array)
%   filename  - Output HDF5 filename (string)
%   celltype  - Cell type information (string array)
%   s         - Additional numeric metadata
%
% Example:
% X = sparse([1 0; 0 2]);
% genelist = ["Gene1", "Gene2"];
% e_writeh5(X, genelist, 'output.h5');

if nargin < 5, batchid = []; end
if nargin < 4, celltype = []; end

if ~ischar(filename) || isempty(filename)
    error('Filename must be a non-empty string.');
end
if ~ismatrix(X)
    error('X must be a 2D matrix.');
end

if nnz(X) == 0
    warning('X is an empty sparse matrix. The output file may not contain meaningful data.');
end


% if isa(X,'SingleCellExperiment')
%     genelist=X.g;
%     X=X.X;
% end
% if nargin<2, genelist=string([1:size(X,1)].'); end

% X = sparse([1 0 2; 0 0 3; 4 5 6]);
% pkg.e_writeh5(X,["a"],'test.h5');
% [Y,g]=pkg.e_readh5('test.h5');


% https://www.10xgenomics.com/support/software/space-ranger/advanced/hdf5-feature-barcode-matrix-format
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html

try
    if issparse(X)
        disp('writing sparse...');
        [indptr, indices, data] = convert_sparse_to_indptr(X);
        h5create(filename, '/shape', size(size(X)));
        h5write(filename, '/shape', size(X));
        
        h5create(filename, '/data', size(data));
        h5write(filename, '/data', data);
    
        h5create(filename, '/indptr', size(indptr));
        h5write(filename, '/indptr', indptr);
        
        h5create(filename, '/indices', size(indices));
        h5write(filename, '/indices', indices);
    else
    
        h5create(filename, '/X', size(X));
        h5write(filename, '/X', X);
    end       
catch ME
        error('Failed to write matrix X to HDF5 file: %s', ME.message);
end

[n, m] = size(X);

if isempty(genelist), genelist = "G"+string(1:n); end
if isempty(celltype), celltype = repmat("undetermined", m, 1); end
if isempty(batchid), batchid = string(ones(m,1)); end

    if ~isempty(genelist)
        h5create(filename, '/g', size(genelist), 'Datatype', 'string');
        h5write(filename, '/g', genelist);
    end
    
    if ~isempty(celltype)
        h5create(filename, '/celltype', size(celltype), 'Datatype', 'string');
        h5write(filename, '/celltype', celltype);
    end
    
    if ~isempty(batchid)
        if ~isstring(batchid)
            batchid = string(batchid);
        end
        h5create(filename, '/batchid', size(batchid), 'Datatype', 'string');
        h5write(filename, '/batchid', batchid);
    end


end



function [indptr, indices, data] = convert_sparse_to_indptr(X)

    % Check if X is sparse
    if ~issparse(X)
        error('Input matrix X must be a sparse matrix');
    end
    
    % Get matrix dimensions
    [~, n] = size(X);
    
    % Initialize indptr with 1 and n+1
    indptr = [1, n+1];
    
    % Find non-zero elements and their indices
    [row, col] = find(X);
    
    % Sort by columns for efficient construction
    [~, sort_idx] = sort(col);
    row = row(sort_idx);
    col = col(sort_idx);
    
    % Accumulate column counts for indptr
    for i = 1:n
        indptr(i+1) = indptr(i) + sum(col == i);
    end
    
    % Assign indices and data
    indices = row;
    % data = full(X(row, col));  % Extract non-zero values
    data = nonzeros(X);

end
