function [idx] = sc_similarg(X, genelist, targetgene, k, isgene)
%KNN network of genes for identifying similar genes
if nargin < 5, isgene = true; end
if nargin < 4, k = 5; end
if nargin < 3
    error('Need 3 inputs (X,genelist,targetgene)');
end

if isnumeric(targetgene)
    targetidx = targetgene;
else
    [y, targetidx] = ismember(targetgene, genelist);
    if ~y, idx = 0;
        return;
    end
end
if ~(targetidx > 0), idx = 0;
    return;
end
if isgene
    u = mean(X, 2, 'omitnan');
    cv2 = var(X, 0, 2, 'omitnan') ./ u.^2;
    A = [u, cv2];
    idx = knnsearch(A, A(targetidx, :), 'K', k+1, ...
        'Distance', 'Euclidean');
    idx = idx(2:end);
end
end
