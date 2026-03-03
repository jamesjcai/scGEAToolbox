function [X, keptidx] = sc_filterc(X, cutoff, verbose)
arguments
    X {mustBeNumeric}
    cutoff (1,1) double {mustBeNonnegative} = 1
    verbose (1,1) logical = false
end

[s] = sum(X, 1);
keptidx = s >= cutoff;
X = X(:, keptidx);

if verbose, fprintf('%d samples (cells) removed.\n', sum(~keptidx)); end

end
