function [ Y ] = mio_mat2cell( X )

[r c]=size(X);

if r>0 & c>0
    Y=mat2cell(X,ones(1,r),ones(1,c));
else
    Y={[]};
end


end

