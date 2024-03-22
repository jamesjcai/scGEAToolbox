function [lgu, dropr, lgcv, genes, X, ...
    removedgeneidx, removedT] = sc_genestat(X, genelist, ...
    sortit, removeinf)

if nargin < 4, removeinf = true; end
if nargin < 3, sortit = true; end
if nargin < 2 || isempty(genelist)
    genelist = "Gene"+string(1:size(X,1)).'; 
end
genelist=genelist(:);

geneidx = 1:length(genelist);
dropr = 1 - sum(X > 0, 2) ./ size(X, 2);

% m = X./sum(X);
% lgm = log(std(m, [], 2, 'omitnan') ./ mean(m, 2, 'omitnan'));


u = mean(X, 2, 'omitnan');
cv = std(X, [], 2, 'omitnan') ./ u;
lgu = log(u+1);
lgcv = log(cv+1);
genes = genelist;


if sortit
    [xyz, si] = sortrows([lgu, dropr, lgcv], [1, 3, 2]);
    lgu = xyz(:, 1);
    dropr = xyz(:, 2);
    lgcv = xyz(:, 3);
    genes = genelist(si); 
    X=X(si,:);
    geneidx = geneidx(si);
end

si = isnan(lgu) | isinf(lgu) | isnan(lgcv) | isinf(lgcv);

if removeinf && any(si)
    removedT = table(genes(si), lgu(si), lgcv(si), dropr(si), ...
        zeros(sum(si),1), ...
        ones(sum(si),1), ...
        ones(sum(si),1));
    lgu(si) = [];
    lgcv(si) = [];
    dropr(si) = [];
    genes(si) = [];
    X(si, :) = [];
    removedgeneidx = geneidx(si);    
else
    removedgeneidx = [];
    removedT = [];
end

end
