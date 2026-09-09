function [T, Xsorted_completed, gsorted_completed, ...
    xyz1] = sc_splinefit(X, genelist, sortit, plotit, removenan)
% SC_SPLINEFIT identify genes with a profile deviated from normal
%
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [X]=sc_norm(X,'type','libsize');
% >> [T]=sc_splinefit(X,genelist,true,true);

% if nargin<2 || isempty(genelist)
%     genelist=string(1:size(X,1))';
%     genelist=strcat("gene_",genelist);
% end

if nargin < 5, removenan = false; end
if nargin < 4, plotit = false; end
if nargin < 3, sortit = true; end
if nargin < 2 || isempty(genelist)
    genelist = string(1:size(X, 1));
end

idx = sum(X, 2) == 0;
if any(idx)
   genelist(idx) = [];
   X(idx, :) = [];
   warning('Empty genes are removed.');
end

m = size(X, 2);
if m < 10000
    [lgu, dropr, lgcv, gsorted, Xsorted, ...
        removedgidx, removedT] = sc_genestat(X, genelist, sortit, removenan);
else
    [lgu, dropr, lgcv, gsorted, Xsorted, ...
        removedgidx, removedT] = pkg.sc_genestat_sparse_blocked(X, genelist, sortit, removenan);
end

if removenan && ~isempty(removedgidx)
    gsorted_completed = [gsorted; genelist(removedgidx)];
    Xsorted_completed = [Xsorted; X(removedgidx, :)];
    assert(isequal(size(Xsorted_completed), size(X)), 'SC_SPLINEFIT')
    assert(length(gsorted_completed) == length(genelist), 'SC_SPLINEFIT')
else
    gsorted_completed = gsorted;
    Xsorted_completed = Xsorted;
end

% lgu=zscore(lgu);
% dropr=zscore(dropr);
% lgcv=zscore(lgcv);

% [~,i]=max(lgcv);

xyz = [lgu, lgcv, dropr];

% [~,j]=sort(pdist2(xyz,xyz(i,:)));
% xyz=xyz(j,:)';
% lgu=lgu(j);
% dropr=dropr(j);
% lgcv=lgcv(j);

% xyz=[lgu dropr lgcv]';

% s = cumsum([0;sqrt(diff(lgu(:)).^2 + diff(dropr(:)).^2 ...
%     + diff(lgcv(:)).^2)]);
s = cumsum([0; sqrt(diff(lgu(:)).^2 + diff(lgcv(:)).^2 ...
    + diff(dropr(:)).^2)]);

% ONCLEANUP rather than a bare off/on pair. The pair left the warning
% disabled for the rest of the session if SPLINEFIT or PPVAL threw between
% the two lines, and its "on" re-enabled a warning the caller may have
% turned off deliberately rather than restoring what they had.
warnState = warning('off', 'MATLAB:rankDeficientMatrix');
restoreWarn = onCleanup(@() warning(warnState));

pp1 = splinefit(s, xyz.', 15, 0.75);
xyz1 = ppval(pp1, s)';

[nearidx, d] = dsearchn(xyz1, xyz);

fitmeanv = xyz1(:, 1);
x = xyz(:, 1); y = xyz(:, 2);
aboveRange = x > max(fitmeanv);
belowRange = x < min(fitmeanv);
belowCurve = (y - xyz1(:, 2)) < 0;
d(aboveRange) = d(aboveRange) ./ 100;
d(belowRange) = d(belowRange) ./ 10;
d(belowCurve) = d(belowCurve) ./ 100;

% D = pdist2(xyz, xyz1);
% d = min(D, [], 2);

% The three lines above are a heuristic down-weighting, not a
% transformation: a gene whose mean falls outside the fitted range, or
% which sits below the curve in CV, is declared not a candidate by dividing
% its distance by 100 or 10. On real 10x data that is 94.8% of genes, and
% their deflated distances sit ~50x below the rest, so any null fitted to
% all of D is fitted to the deflation rather than to the noise.
%
% So the null is fitted on the candidates, and everything else gets p = 1.
% An arbitrary constant times a distance has no null, and reporting one for
% it was most of what went wrong here: the old code fitted a symmetric
% null to the whole of D by mirroring it, then read a single tail of that
% null with a scale taken from the standard deviation of a truncated bulk.
% Every one of those choices pushed the same way. On the bundled 10x
% example it called 23.3% of all genes significant at p<0.05.
%
% The D ranking is untouched, so gene selection by D -- which is what every
% caller of this function actually uses; nothing reads PVAL or FDR -- is
% exactly as before.
isCandidate = ~(aboveRange | belowRange | belowCurve);
pval = ones(size(d));
if any(isCandidate)
    pval(isCandidate) = pkg.e_deviationpvalue(d(isCandidate));
end
fdr = pkg.e_fdr(pval);

if ~isempty(gsorted)
    genes = gsorted;
    T = table(genes, lgu, lgcv, dropr, d, pval, fdr, nearidx);
else
    T = table(lgu, lgcv, dropr, d, pval, fdr, nearidx);
end
% 'variablenames',{'Genes','Log10_Mean','Dropout_Rate','Log10_CV','Deviation_3DFeature'});

% T.d(T.dropr > (1 - 0.05)) = 0; % ignore genes with dropout rate > 0.95
% T.d(T.dropr < (0.01)) = 0;     % ignore genes with dropout rate < 0.01 (removes ribosomal and mitochondrial genes)

% disp('NOTE: Genes with dropout rate > 0.95 are excluded.');

if ~isempty(removedT) && istable(removedT)
    removedT.Properties.VariableNames = T.Properties.VariableNames;
    T = [T; removedT];
end

if sortit
    [T, idx] = sortrows(T, 'd', 'descend');
    gsorted_completed = gsorted_completed(idx);
    Xsorted_completed = Xsorted_completed(idx, :);
end

if length(gsorted_completed) ~= length(genelist)
    error('Output GENES are less than input GENES (some GENES are removed).');
end

if plotit
    figure;
    scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'filled', 'MarkerFaceAlpha', .1);
    hold on
    plot3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4);
    xlabel('Mean, log');
    ylabel('CV, log');
    zlabel('Dropout rate (% of zeros)');

    if ~isempty(gsorted)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn3, gsorted, Xsorted};
    end
    hold off
end

end
