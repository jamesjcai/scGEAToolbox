function [X, genelist] = e_readh5(filenm)

% https://www.10xgenomics.com/support/software/space-ranger/advanced/hdf5-feature-barcode-matrix-format
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html

% X = sparse([1 0 2; 0 0 3; 4 5 6]);
% pkg.e_writeh5(X,["a"],'test.h5');
% [Y,g]=pkg.e_readh5('test.h5');

    try
        X = h5read(filenm, '/X');
    catch
        disp('reading sparse...');
        grouptag = '/';
        data = pkg.e_guessh5field(filenm, {grouptag}, {'data'}, true);
        indices = pkg.e_guessh5field(filenm, {grouptag}, {'indices'}, true);
        indptr = pkg.e_guessh5field(filenm, {grouptag}, {'indptr'}, true);
        shape = pkg.e_guessh5field(filenm, {grouptag}, {'shape'}, true);

        X = spalloc(shape(1), shape(2), length(data));
        for k = 1:length(indptr)-1
            ix = (indptr(k) + 1:indptr(k+1))-1;
            X(indices(ix), k) = data(ix);
        end
    end
    genelist = h5read(filenm, '/g');
end


% https://stackoverflow.com/questions/43021896/construct-sparse-matrix-in-matlab-from-compressed-sparse-column-csc-format