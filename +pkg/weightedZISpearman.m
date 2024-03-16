%' the weightedZISpearman function calculates weighted rho\*,
%' where rho\* is described in Pimentel et al (2009).
%' This association measure is defined for zero-inflated,
%' non-negative random variables.
%'
%' @title weightedZISpearman
%' @param w weight vector, values should be between 0 and 1
%' @param x x and y are non-negative data vectors
%' @param y x and y are non-negative data vectors
%' @return \code{numeric} weighted rho* association value between x and y
%'
%'  Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho for
%'  Zero-Inflated Data" (2009). Dissertations. 721.
%'  https://scholarworks.wmich.edu/dissertations/721

%' @examples
%'
%' x = pmax(0,rnorm(100))
%' y = pmax(0,rnorm(100))
%' w = runif(100)
%' weightedZISpearman(x,y,w)
%'
%' @export
%{
cd E:\RStudio
x=readmatrix('xx.csv'); x=x(:,2);
y=readmatrix('yy.csv'); y=y(:,2);
w=readmatrix('ww.csv'); w=w(:,2);
r=-0.07823878;
%}

function [rho_star] = weightedZISpearman(x, y, w)
if nargin < 3, w = 1.0; end
if any(x < 0 | y < 0), error('x and/or y values have negative values'); end
if (length(x) ~= length(y)), error('x and y should have the same length'); end
if (isscalar(w)), w = w * ones(size(x)); end

i = x > 0;
j = y > 0;
k = i & j;

p_11 = sum(w.*k) / sum(w);
p_00 = sum(w.*(~i & ~j)) / sum(w);
p_01 = sum(w.*(~i & j)) / sum(w);
p_10 = sum(w.*(i & ~j)) / sum(w);
rho_11 = pkg.weightedSpearman(x(k), y(k), w(k));
rho_star = p_11 * (p_01 + p_11) * (p_10 + p_11) * rho_11 + ...
    3 * (p_00 * p_11 - p_10 * p_01);
if isnan(rho_star)
    % warning("Zero inflated Spearman correlation is undefined, returning Spearman correlation");
    rho_star = pkg.weightedSpearman(x, y, w);
end
end
