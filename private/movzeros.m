function [y] = movzeros(x, k)
if nargin < 2, k = 20; end
y = myFilter(x, k, @(x)sum(x == 0));
end


function y = myFilter(x, k, fcn)
% https://code.i-harness.com/en/q/79aed8
idx = hankel(1:k, k:length(x));
y = cellfun(fcn, num2cell(x(idx), 1));
end
