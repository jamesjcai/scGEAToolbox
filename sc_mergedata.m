function [X, genelist, c] = sc_mergedata(X1, X2, genelist1, genelist2, method, ignorecase)
if nargin < 6, ignorecase = false; end
if nargin < 5, method = 'intersect'; end
if ignorecase
    genelist1 = upper(genelist1);
    genelist2 = upper(genelist2);
end

switch lower(method)
    case 'intersect'
        [genelist, i, j] = intersect(genelist1, genelist2, 'stable');
        X1a = X1(i, :);
        X2a = X2(j, :);
    case 'union'
        [genelist] = union(genelist1, genelist2, 'stable');
        [~, i] = ismember(genelist1, genelist);
        [~, j] = ismember(genelist2, genelist);
        if issparse(X1) || issparse(X2)
            X1a = sparse(length(genelist), size(X1, 2));
            X2a = sparse(length(genelist), size(X2, 2));
        else
            X1a = zeros(length(genelist), size(X1, 2));
            X2a = zeros(length(genelist), size(X2, 2));
        end
        X1a(i, :) = X1;
        X2a(j, :) = X2;
end
X = [X1a, X2a];
c = 1 + [zeros(size(X1, 2), 1); ones(size(X2, 2), 1)];
end