function [T, sx, sy, sz, d, genesfit] = sc_splinefit2(X, Y, genelistx, genelisty, sortid)
%SC_SPLINEFIT2  Compare two datasets' gene profiles by spline fit.
%
%   T = SC_SPLINEFIT2(X, Y, genelistx, genelisty) fits SC_SPLINEFIT to each
%   dataset over their shared genes and joins the two tables, adding
%   T.dd = d2 - d1, the change in each gene's deviation from its own
%   dataset's fitted curve.
%
%   [T, sx, sy, sz, d, genesfit] = SC_SPLINEFIT2(...) also returns the two
%   fitted curves, the Procrustes transform of SY onto SX, and the
%   Procrustes dissimilarity. SX, SY and SZ share one gene order, returned
%   as GENESFIT.
%
%   SORTID (default true) sorts T by dd descending.

if nargin < 3, error('X,Y,GENELIST are required.'); end
if nargin < 4, genelisty = genelistx; end
if nargin < 5, sortid = true; end

[genelist, i, j] = intersect(genelistx, genelisty, 'stable');
X = X(i, :);
Y = Y(j, :);

% The third output is the gene order each fit came back in, and it is
% needed: SC_SPLINEFIT sorts genes through SC_GENESTAT by that dataset's
% own statistics, sortrows([lgu, dropr, lgcv], [1, 3, 2]), so row r of SX
% and row r of SY are different genes. The INTERSECT above exists precisely
% to make row r the same gene in X and Y, and the sort inside undoes it.
[T1, ~, g1, sx] = sc_splinefit(X, genelist);

T1.Properties.VariableNames = {'genes', 'logu1', 'logcv1',...
    'dropr1', 'd1', 'pval1', 'fdr1', 'nearidx1'};
[T2, ~, g2, sy] = sc_splinefit(Y, genelist);
T2.Properties.VariableNames = {'genes', 'logu2', 'logcv2',...
    'dropr2', 'd2', 'pval2', 'fdr2', 'nearidx2'};

T = join(T1, T2, 'Keys', 'genes');
T.dd = T.d2 - T.d1;
if sortid
    T = sortrows(T, 'dd', 'descend');
end

% Re-pair the two curves by gene before comparing them. This used to be
% procrustes(sx, sy) straight off the two fits, matching them by row index.
% Measured on 300 genes in two conditions: only 2 of the 300 rows named the
% same gene in both orders, and the dissimilarity came out at 0.012307
% against 0.995189 once paired by name. The error is in the flattering
% direction, which is why it could sit unnoticed: both curves are the same
% spline shape traced through the same statistic space, so Procrustes fits
% one to the other very well whatever gene sits in which row, and two
% genuinely different datasets were reported as nearly identical.
%
% INTERSECT rather than an assert, because SC_SPLINEFIT also drops genes
% that are empty in its own input, so the two fits can have different gene
% sets: one gene empty in X alone made procrustes fail outright with
% stats:procrustes:InputSizeMismatch. SC_DVG re-establishes the same
% correspondence for the same reason, by sorting both tables on gene name.
[genesfit, ia, ib] = intersect(string(g1), string(g2), 'stable');
sx = sx(ia, :);
sy = sy(ib, :);

[d, sz] = procrustes(sx, sy);
