function y = nanmedian(x,dim)
%NANMEDIAN Median value, ignoring NaNs.
%   M = NANMEDIAN(X) returns the sample median of X, treating NaNs as
%   missing values.  For vector input, M is the median value of the non-NaN
%   elements in X.  For matrix input, M is a row vector containing the
%   median value of non-NaN elements in each column.  For N-D arrays,
%   NANMEDIAN operates along the first non-singleton dimension.
%
%   NANMEDIAN(X,'all') is the median value of all the elements of X.
%
%   NANMEDIAN(X,DIM) takes the median along the dimension DIM of X.
%
%   NANMEDIAN(X,VECDIM) finds the median values of the elements of X based 
%   on the dimensions specified in the vector VECDIM.
%
%   See also MEDIAN, NANMEAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2018 The MathWorks, Inc.


if nargin == 1
    y = prctile(x, 50);
else
    y = prctile(x, 50,dim);
end
