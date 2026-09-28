function [X, g, b, c, filenm] = sc_read10xh5file(filenm)
% Read 10x Genomics H5 file
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices

% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3489183
% h5file='GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5';

X = []; % expression matrix
g = []; % gene names
b = []; % barcode of cells
c = []; % batch id
if nargin < 1 || isempty(filenm)
    [filenm, pathname] = uigetfile( ...
        {'*.h5;*.hdf5', 'HDF5 Files (*.h5)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a 10x Genomics H5 file');
    if isequal(filenm, 0), return; end
    filenm = fullfile(pathname, filenm);
end
if exist(filenm, 'file') ~= 2
    error('FileNotFound');
end

grouptag = "/matrix/";
    h = h5info(filenm);
    assert(any(contains(string({h.Groups.Name}), "/matrix")))

data = pkg.e_guessh5field(filenm, {grouptag}, {'data'}, true);
indices = pkg.e_guessh5field(filenm, {grouptag}, {'indices'}, true);
indptr = pkg.e_guessh5field(filenm, {grouptag}, {'indptr'}, true);
shape = pkg.e_guessh5field(filenm, {grouptag}, {'shape'}, true);

g = pkg.e_guessh5field(filenm, {grouptag, '/matrix/features/'}, {'gene_names', 'name'}, false);
if isempty(g)
    warning('Gene names or feature names are not assigned.');
end

b = pkg.e_guessh5field(filenm, {grouptag, '/matrix/features/'}, {'barcodes'}, false);
if isempty(b), warning('B is not assigned.'); end

if ~isMATLABReleaseOlderThan('R2025a')
    X = spalloc(shape(1), shape(2), length(data), 'single');
else
    X = spalloc(shape(1), shape(2), length(data));
end

for k = 1:length(indptr) - 1
    ix = indptr(k) + 1:indptr(k+1);
    X((indices(ix) + 1), k) = data(ix);
end

g = deblank(string(g));

if all(contains(b,'-'))
    c = extractAfter(b, "-");
end

end

%{

% function countMatrix = getMatrixFromH5(filename)
%     info = h5info(filename, '/matrix');
%
%     barcodes = h5read(filename, '/matrix/barcodes');
%     data = h5read(filename, '/matrix/data');
%     indices = h5read(filename, '/matrix/indices');
%     indptr = h5read(filename, '/matrix/indptr');
%     shape = h5read(filename, '/matrix/shape');
%
%     matrix = sparse(indices+1, indptr+1, data, shape(2), shape(1));
%
%     featureRef = struct();
%     featureGroup = info.Groups(strcmp({info.Groups.Name}, '/matrix/features'));
%     featureRef.id = h5read(filename, '/matrix/features/id');
%     featureRef.name = h5read(filename, '/matrix/features/name');
%     featureRef.featureType = h5read(filename, '/matrix/features/feature_type');
%
%     tagKeys = h5read(filename, '/matrix/features/_all_tag_keys');
%     for i = 1:length(tagKeys)
%         key = char(tagKeys(i));
%         featureRef.(key) = h5read(filename, ['/matrix/features/' key]);
%     end
%
%     countMatrix = struct('featureRef', featureRef, 'barcodes', barcodes, 'matrix', matrix);
% end
%
% filteredH5 = 'filtered_feature_bc_matrix.h5';
% filteredMatrix = getMatrixFromH5(filteredH5);

% https://www.10xgenomics.com/support/software/space-ranger/advanced/hdf5-feature-barcode-matrix-format

%}
