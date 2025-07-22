function M=divide_sparse(M, divisor)
%DIVIDE_SPARSE Divide each row in the sparse matrix M by the number in the
% corresponding row of "divisor". In practice, this is much faster for large
% arrays than using normal MATLAB syntax.
%
% M = DIVIDE_SPARSE(M, divisor)
%
% Parameters
% ----------
% M: sparse matrix of size (m1, m2)
% 
% divisor: array of size (m1, 1)
% 
% Returns
% -------
% M: sparse matrix of size (m1, m2)
%     The result of dividing the (i, j)-th entry of M by the i-th component of
%     "divisor".

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

    R=size(M, 1);
    divisor = full(divisor);
    
    M = spdiags(1 ./ divisor, 0, R, R) * M;
end