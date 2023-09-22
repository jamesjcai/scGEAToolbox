% https://github.com/shazanfar/DCARS

%' the weightedSpearman function
%'
%' @title weightedSpearman
%' @param x x and y are data vectors
%' @param y x and y are data vectors
%' @param w weight vector, values should be between 0 and 1
%' @return \code{numeric} weighted correlation value between x and y

%' @examples
%'
%' x = rnorm(100)
%' y = rnorm(100)
%' w = runif(100)
%' weightedSpearman(x,y,w)
%'
%' @export

function [r] = weightedSpearman(x, y, w)
if nargin < 3, w = 1.0; end
if any(x < 0 | y < 0), error('x and/or y values have negative values'); end
if (length(x) ~= length(y)), error('x and y should have the same length'); end
if (length(w) == 1), w = w * ones(size(x)); end

i = w > 0;
% [~, ~, rank] = unique(Data);
xr = tiedrank(x(i));
yr = tiedrank(y(i));
r = pkg.weightedPearson(xr, yr, w(i));
end