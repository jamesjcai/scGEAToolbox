function npary = mat2np(mat)

% convert matlab matrix to python (Numpy) ndarray 
sh = fliplr(size(mat));
mat2 = reshape(mat,1,numel(mat));  % [1, n] vector
npary = py.numpy.array(mat2);
npary = npary.reshape(int32(sh)).transpose();  % python ndarray

end

% https://gist.github.com/samtx
