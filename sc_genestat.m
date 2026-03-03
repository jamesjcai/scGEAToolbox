function [lgu, dropr, lgcv, genes, X, removedIdx, removedT] = sc_genestat(X, genelist, sortit, removeinf)
% SC_GENESTAT  Compute per-gene statistics and optionally filter invalid values
%
%   [lgu,dropr,lgcv,genes,X] = sc_genestat(X)
%   [lgu,dropr,lgcv,genes,X] = sc_genestat(X,genelist,sortit,removeinf)
%
%   Inputs:
%     X         m-by-n gene expression matrix (genes × cells)
%     genelist  m×1 string or cellstr vector of gene names
%     sortit    true|false – sort output by (lgu, lgcv, dropr). Default true
%     removeinf true|false – remove NaN/Inf entries. Default true
%
%   Outputs:
%     lgu        log1p(mean expression)
%     dropr      dropout rate (fraction of zeros)
%     lgcv       log1p(coefficient of variation)
%     genes      gene names (possibly sorted/filtered)
%     X          filtered (and re-ordered) expression matrix
%     removedIdx indices of removed genes (if any)
%     removedT   table of removed gene stats (if any)

arguments
    X {mustBeNumeric, mustBeNonempty}
    genelist = []
    sortit (1,1) logical = true
    removeinf (1,1) logical = true
end

% --- Setup gene names ---
nGenes = size(X, 1);
if isempty(genelist)
    genelist = pkg.i_defaultgenenames(nGenes);
else
    if iscellstr(genelist)
        genelist = string(genelist);
    end
    genelist = genelist(:);
end

% --- Compute stats ---
u = mean(X, 2, 'omitnan');
cv = std(X, 0, 2, 'omitnan') ./ u;
lgu = log1p(u);
lgcv = log1p(cv);
dropr = 1 - sum(X > 0, 2) ./ size(X, 2);
genes = genelist;
idx = (1:nGenes).';

% --- Sort rows if requested ---
if sortit
    M = [lgu, dropr, lgcv];
    [~, order] = sortrows(M, [1, 3, 2]);
    lgu = lgu(order);
    dropr = dropr(order);
    lgcv = lgcv(order);
    genes = genes(order);
    X = X(order, :);
    idx = idx(order);
end

% --- Identify and optionally remove invalid entries ---
mask = isnan(lgu) | isinf(lgu) | isnan(lgcv) | isinf(lgcv);
removedIdx = [];
removedT = table();

if removeinf && any(mask)
    removedT = table(genes(mask), lgu(mask), lgcv(mask), dropr(mask), ...
                     'VariableNames', {'Gene','lgu','lgcv','dropout'});
    removedIdx = idx(mask);
    lgu(mask) = [];
    dropr(mask) = [];
    lgcv(mask) = [];
    genes(mask) = [];
    X(mask, :) = [];
end
end