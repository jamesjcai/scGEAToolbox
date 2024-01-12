function [lgu, dropr, lgcv, genes, X, removedgeneidx, removedT] = sc_genestat(X, genelist, sortit, removeinf)

if nargin < 4, removeinf = true; end
if nargin < 3, sortit = true; end
if nargin < 2 || isempty(genelist), genelist = "Gene"+string(1:size(X,1)).'; end

geneidx = 1:length(genelist);

dropr = 1 - sum(X > 0, 2) ./ size(X, 2);
u = mean(X, 2, 'omitnan');
cv = std(X, [], 2, 'omitnan') ./ u;
lgu = log(u);
lgcv = log(cv);
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

if removeinf
    si = isnan(lgu) | isinf(lgu) | isnan(lgcv) | isinf(lgcv);

    removedT = table(genes(si), lgu(si), lgcv(si), dropr(si), zeros(size(dropr(si))), ...
            ones(size(dropr(si))), ones(size(dropr(si))) );


    lgu(si) = [];
    lgcv(si) = [];
    dropr(si) = [];
    genes(si) = [];
    X(si, :) = [];
    if length(genes) ~= length(genelist)
        warning('Output GENES are less than input GENES (some GENES are removed).');
    end
    removedgeneidx = geneidx(si);
    
else
    removedgeneidx = [];
    removedT = [];
end

end
