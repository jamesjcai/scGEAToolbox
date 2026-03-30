function [T, X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2] = ...
        sc_dvg(sce1, sce2, cL1, cL2, method)
% SC_DVG - Differential variability (DV) analysis between two groups
%
% Inputs:
%   sce1, sce2 : SingleCellExperiment objects for each group
%   cL1, cL2   : group label cell arrays (optional, default {'1'},{'2'})
%   method     : 'splinefit' (default, PMID:40113778) | 'brennecke' (PMID:24056876)
%
% Outputs:
%   T    : results table sorted by DiffDist descending
%   X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2 :
%          splinefit visualization data (empty for 'brennecke')

if nargin < 3 || isempty(cL1), cL1 = {'1'}; end
if nargin < 4 || isempty(cL2), cL2 = {'2'}; end
if nargin < 5 || isempty(method), method = 'splinefit'; end

X1 = []; X2 = []; g = []; xyz1 = []; xyz2 = [];
px1 = []; py1 = []; pz1 = []; px2 = []; py2 = []; pz2 = [];

if sce1.NumCells < 50 || sce2.NumCells < 50
    warning('One of groups contains too few cells (n < 50). The result may not be reliable.');
end
if sce1.NumGenes < 50 || sce2.NumGenes < 50
    warning('One of groups contains too few genes (n < 50). The result may not be reliable.');
end

if ~isequal(sce1.g, sce2.g)
    [g_ori, ia, ib] = intersect(sce1.g, sce2.g, 'stable');
    X1_ori = sce1.X(ia, :);
    X2_ori = sce2.X(ib, :);
else
    g_ori = sce1.g;
    X1_ori = sce1.X;
    X2_ori = sce2.X;
end

switch lower(method)
    case 'splinefit'
        X1_ori = sc_norm(X1_ori, 'type', 'libsize');
        X2_ori = sc_norm(X2_ori, 'type', 'libsize');

        [T1, X1, g1, xyz1] = sc_splinefit(X1_ori, g_ori, true, false);
        [T1, idx1] = sortrows(T1, 'genes', 'ascend');
        X1 = X1(idx1, :);
        g1 = g1(idx1);

        [T2, X2, g2, xyz2] = sc_splinefit(X2_ori, g_ori, true, false);
        [T2, idx2] = sortrows(T2, 'genes', 'ascend');
        X2 = X2(idx2, :);
        g2 = g2(idx2);

        assert(isequal(g1, g2));
        g = g1;

        px1 = T1.lgu; py1 = T1.lgcv; pz1 = T1.dropr;
        px2 = T2.lgu; py2 = T2.lgcv; pz2 = T2.dropr;

        v1 = ([px1 py1 pz1] - xyz1(T1.nearidx, :));
        v2 = ([px2 py2 pz2] - xyz2(T2.nearidx, :));

        DiffDist = vecnorm(v1 - v2, 2, 2);
        DiffSign = sign(vecnorm(v1, 2, 2) - vecnorm(v2, 2, 2));

        ddz = zscore(DiffDist);
        pval = 1 - normcdf(ddz);
        clear ddz;

        T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s', cL1{1}));
        T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s', cL2{1}));
        idxx = T1.(8) == 1 | T2.(8) == 1 | T1.(8) == max(T1.(8)) | T2.(8) == max(T2.(8));
        T1 = T1(:, 1:end-1);
        T2 = T2(:, 2:end);
        T2 = T2(:, 1:end-1);
        T1.Properties.VariableNames{1} = 'gene';

        T = [T1 T2 table(DiffDist) table(DiffSign) table(pval)];
        T.DiffDist(idxx) = 0;
        T = sortrows(T, 'DiffDist', 'descend');

    case 'brennecke'
        [T1] = sc_hvg(X1_ori, g_ori, true, false);
        T1(T1.removedidx1 | T1.removedidx2, :) = [];
        [T2] = sc_hvg(X2_ori, g_ori, true, false);
        T2(T2.removedidx1 | T2.removedidx2, :) = [];
        [gene, idx1, idx2] = intersect(T1.genes, T2.genes);

        DiffDist = T1.residualcv2(idx1) - T2.residualcv2(idx2);
        DiffDistAbs = abs(DiffDist);
        DiffSign = sign(DiffDist);
        ddz = zscore(DiffDistAbs);
        pval = 1 - normcdf(ddz);

        T = table(gene, DiffDist, DiffDistAbs, DiffSign, pval);
        T = sortrows(T, 'DiffDistAbs', 'descend');

    otherwise
        error('sc_dvg: unknown method ''%s''. Use ''splinefit'' or ''brennecke''.', method);
end
end
