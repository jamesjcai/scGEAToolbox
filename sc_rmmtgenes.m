function [X, genelist, idx] = sc_rmmtgenes(X, genelist, mtgenenamepat, verbose)
%Remove mt-genes
if nargin < 3
    mtgenenamepat = "mt-";
end
if nargin < 4, verbose = true; end
idx = startsWith(genelist, mtgenenamepat, 'IgnoreCase', true);
if any(idx)
    if verbose
        genelist(idx);
        fprintf('%d mt-genes found and removed.\n', sum(idx));
    end
    genelist = genelist(~idx);
    X = X(~idx, :);
else
    if verbose
        fprintf('No mt-genes found or removed.\n');
    end
end
if nargout > 3
    idx = find(idx);
end

% find(contains(genelist,"ENSG00000198804"))  % MT-CO1
end