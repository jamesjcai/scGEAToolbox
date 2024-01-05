function [X, genelist] = e_readh5(filenm)


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
        size(X)
        for k = 1:length(indptr)-1
            ix = indptr(k) + 1:indptr(k+1);
            ix = ix - 1;
            y = indices(ix);
            X(y, k) = data(ix);
        end
    end
    genelist = h5read(filenm, '/g');
end
