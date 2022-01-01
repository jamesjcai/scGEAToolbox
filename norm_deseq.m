function [x,sf]=norm_deseq(x)
%DESeq normalization
% For DESeq normalization, the geometric mean for each gene was computed 
% after removing all zeroes. This is necessary to avoid a situation where 
% a majority of genes have geometric means of zero, such that the majority
% of ratios to the geometric mean would be undefined. Size factors were 
% then computed using the estimateSizeFactorsForMatrix function in DESeq2 
% v1.10.1 [21]. In this function, ratios of zero were automatically removed
% prior to calculation of the median in each library, to avoid obtaining a
% size factor equal to zero        
    y=nan(size(x));
    y(x>0)=x(x>0);
    sf=median(y./nangeomean(y,2),'omitnan');
    x=x./sf;
end


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
end