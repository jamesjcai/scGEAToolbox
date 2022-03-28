% R = matthewscorr(X, Y)
%
% computes the 'Matthews correlation coefficient' as described on Wikipedia
% (http://en.wikipedia.org/wiki/Matthews_correlation_coefficient)
%
% It's supposedly a special correlation coefficient for binary data (e.g.
% in binary classification) and computes the correlation based on true and
% false postivies and negatives, but the correlation coefficients seem to
% be the same as those computed by Matlab's corrcoef.
%
% in:
%       X   -   binary data set 1
%               [N, nsets] = size
%       Y   -   binary data set 2
%               [N, nsets] = size
% out:
%       R   -   Matthews correlation coefficient between columns of X and Y
%               [1, nsets] = size
% author:
%       Sebastian Bitzer (bitzer@cbs.mpg.de)
function R = matthewscorr(X, Y)

TP = sum(X == 1 & Y == 1);
FP = sum(X == 0 & Y == 1);
TN = sum(X == 0 & Y == 0);
FN = size(X, 1) - TP - FP - TN;

R = (TP .* TN - FP .* FN) ./ ...
    sqrt( (TP + FP) .* (TP + FN) .* (TN + FP) .* (TN + FN) );