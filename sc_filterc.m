function [X, keptidx] = sc_filterc(X, cutoff, verbose)
if nargin < 2, cutoff = 1; end
if nargin < 3, verbose = false; end
[s] = sum(X, 1);
keptidx = s >= cutoff;
X = X(:, keptidx);
if verbose, fprintf('%d samples (cells) removed.\n', sum(~keptidx)); end
end
