function [X, genelist] = sc_filterg(X, genelist, cutoff, verbose)
arguments
    X {mustBeNumeric}
    genelist = []
    cutoff (1,1) double {mustBePositive} = 1
    verbose (1,1) logical = false
end
% cutoff < 1:  interpreted as minimum non-zero fraction (dropout-rate mode)
% cutoff >= 1: interpreted as minimum total read count (count-threshold mode)

if cutoff < 1.0 % e.g., 0.9
    if verbose
        fprintf('Discard genes with poor expression values (more than %d%% zeros in all cells).\n', ...
            100*cutoff);
    end
    r = sum(X ~= 0, 2) ./ size(X, 2);
    i = r >= cutoff;
else
    [u] = sum(X, 2);
    i = u >= cutoff;
    if verbose
        fprintf('Discard genes with poor expression values (with less than %d reads among all cells).\n', ...
            cutoff);
    end
end
% We discard cells with poor gene expression values (more than 90% zeros in all cells)
% As the default, we filter all the genes with less than 5 reads among 99% of the samples
X = X(i, :);
if ~isempty(genelist), genelist = genelist(i); end
if verbose, fprintf('%d genes removed.\n', sum(~i)); end
end