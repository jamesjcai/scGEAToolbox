function [X, genelist, barcodes, c, filenm] = sc_readhdf5file(filenm)
%Read HDF5 file
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices

% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3489183
% h5file='GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5';
if nargin < 1, filenm = []; end

[X, genelist, barcodes, c, filenm] = sc_read10xh5file(filenm);