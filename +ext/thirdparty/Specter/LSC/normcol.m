% Author: Modified package by Van Hoan Do
% Creates a matrix W from A where the each column is normalized by its sum.
% If A is sparse W is sparse.
%
% Author: Frank Lin (frank@cs.cmu.edu)

function W=normcol(A)

    W=normrow(A')';

end
