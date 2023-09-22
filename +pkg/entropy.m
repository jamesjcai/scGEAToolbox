function [entropy, bins] = entropy(x, k, t)
%ENTROPY - a function to calculate the entropy of a vector using the
%B-spline approximation described in Daub et. al.
%
% Given a real column vector x, spline order k, and a knot t representing
% the binning scheme, calculate the approximate entropy as described in
% [1]. The data x is scaled to the interval [0 (length(t)-2k+1)]. The knots
% described in [1] for estimation are of the form
%
% t = [zeros(1,k-1), 0:(M-k+1) , (M-k+1)*ones(1,k-1)];
%
% where M = length(t)-k, although more general knots could be used with
% this function, provided they have a minimum of 0 and a max of M-k+1.
%
% This is a bioinformatics helper function, and in general, should not be
% called directly by users.
%
% The estimator in [1] is subject to bias. Choice of t and k strongly
% affect this bias in ways discussed in [1], but also, depending on the
% underlying distribution of x and y, the number of samples from x and y
% affect this bias. For the estimate of mutual information to converge in
% a predictable way, the underlying distribution of x and y should be bounded
% above and below.
%
% It's also important to realize that the estimator converges to the mutual
% information of two discrete, proxy random variable for x that is
% defined by the binning process. If x is sampled from a continuous
% random variable, the estimated entropy is not, in general, the
% entropy of the random variable.
%
% References:
% [1] Daub, Carsten O., et al. "Estimating mutual information using
%     B-spline functions?an improved similarity measure for analysing gene
%     expression data." BMC bioinformatics 5.1 (2004): 118.

% Copyright MathWorks 2014


% validate x
if ~(isnumeric(x) && isreal(x) && iscolumn(x))
    error('bioinfo:entropy:InvalidX', 'x must be a real column vector');
end

%% scale x to the domain of the spline

numBins = length(t) - k;

xmin = min(x);
xmax = max(x);

if xmax == xmin
    error('bioinfo:entropy:InvalidXRange', 'x must have more than 1 value')
end

small = eps(xmax-xmin);

z = (x - xmin) * (numBins - k + 1) / (xmax - xmin + small);

z = z(:);

%%

bins = zeros(length(z), numBins);

for iter = 1:numBins
    bins(:, iter) = spval(spmak(t(iter:(iter + k)), 1), z);
end

prob = sum(bins, 1) / length(z); % probability of bin membership for z
entropy = -sum(prob.*log(prob+small));