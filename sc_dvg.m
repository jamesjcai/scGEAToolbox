function [T, X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2] = ...
        sc_dvg(sce1, sce2, cL1, cL2, method, direction)
% SC_DVG - Differential variability (DV) analysis between two groups
%
% Inputs:
%   sce1, sce2 : SingleCellExperiment objects for each group
%   cL1, cL2   : group label cell arrays (optional, default {'1'},{'2'})
%   method     : 'splinefit' (default, PMID:40113778) | 'analytic' |
%                'brennecke' (PMID:24056876)
%
%                'analytic' replaces the fitted smoothing spline with the
%                closed-form gamma-Poisson reference curve of
%                SC_ANALYTICFIT, scoring genes the same way against it. The
%                curve is defined at every mean, so unlike 'splinefit' it
%                has no end-of-fit region and discards no genes there.
%   direction  : how DiffSign, and so up/down, is decided:
%                'mean' (default) - +1 when the gene's mean (library-size
%                    normalized) is higher in sce1 (tested) than in sce2
%                    (baseline): up-regulated.
%                'deviation' - +1 when the gene deviates further from its
%                    curve in sce1 than in sce2: more variable in sce1.
%                    For 'brennecke', the sign of the residual CV^2
%                    difference.
%
% Outputs:
%   T    : results table sorted by DiffDist descending. DiffDist is the
%          size of the variability difference; DiffSign is its direction
%          as chosen by DIRECTION. Split DV genes on DiffSign.
%   X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2 :
%          curve visualization data ('splinefit' and 'analytic'; empty for
%          'brennecke')

if nargin < 3 || isempty(cL1), cL1 = {'1'}; end
if nargin < 4 || isempty(cL2), cL2 = {'2'}; end
if nargin < 5 || isempty(method), method = 'splinefit'; end
if nargin < 6 || isempty(direction), direction = 'mean'; end
direction = validatestring(direction, {'mean', 'deviation'}, ...
    mfilename, 'direction', 6);

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

        % The spline spans only the genes it was fitted through, so genes
        % sitting on its first or last point are dropped -- see IN_SCOREDV.
        [T, px1, py1, pz1, px2, py2, pz2] = ...
            in_scoredv(T1, T2, xyz1, xyz2, cL1, cL2, true, direction);

    case 'analytic'
        % SC_ANALYTICFIT works out alpha and the dropout law from the
        % library sizes, so it takes RAW counts and normalizes internally.
        % X1_ori/X2_ori are therefore passed through unnormalized, and the
        % matrices handed back for plotting are normalized separately so
        % that they match what the 'splinefit' branch returns.
        [T1, xyz1] = sc_analyticfit(X1_ori, g_ori, SortIt=false);
        [T2, xyz2] = sc_analyticfit(X2_ori, g_ori, SortIt=false);

        % SC_ANALYTICFIT returns one foot point per gene, in gene order.
        % Callers draw XYZ as a polyline and expect it to trace the curve,
        % which gene order does not: on the bundled data it steps backwards
        % 4,740 times, by up to 0.51 of a 7.67 range.
        [xyz1, T1.nearidx] = in_ordercurve(xyz1, T1.nearidx);
        [xyz2, T2.nearidx] = in_ordercurve(xyz2, T2.nearidx);

        % Each fit drops the genes that are all-zero in its own sample, so
        % the two tables need not cover the same genes. INTERSECT puts both
        % in one order; NEARIDX is carried along and still indexes into the
        % unsubset XYZ1/XYZ2.
        [~, ia, ib] = intersect(T1.genes, T2.genes, 'stable');
        T1 = T1(ia, :);
        T2 = T2(ib, :);
        [T1, ord] = sortrows(T1, 'genes', 'ascend');
        T2 = T2(ord, :);
        assert(isequal(T1.genes, T2.genes));
        g = T1.genes;

        X1 = sc_norm(X1_ori, 'type', 'libsize');
        X2 = sc_norm(X2_ori, 'type', 'libsize');
        [~, loc] = ismember(g, g_ori);
        X1 = X1(loc, :);
        X2 = X2(loc, :);

        % No end-of-curve exclusion: the analytic curve is defined at every
        % mean, so no gene can fall off the end of it.
        [T, px1, py1, pz1, px2, py2, pz2] = ...
            in_scoredv(T1, T2, xyz1, xyz2, cL1, cL2, false, direction);

    case 'brennecke'
        [T1] = sc_hvg(X1_ori, g_ori, true, false);
        T1(T1.removedidx1 | T1.removedidx2, :) = [];
        [T2] = sc_hvg(X2_ori, g_ori, true, false);
        T2(T2.removedidx1 | T2.removedidx2, :) = [];
        [gene, idx1, idx2] = intersect(T1.genes, T2.genes);

        DiffDist = T1.residualcv2(idx1) - T2.residualcv2(idx2);
        DiffDistAbs = abs(DiffDist);
        % See DIRECTION. SC_HVG normalizes by library size before taking
        % U, so the two means are comparable. DIFFDIST keeps the sign of
        % the residual CV^2 difference either way.
        if strcmp(direction, 'mean')
            DiffSign = sign(T1.u(idx1) - T2.u(idx2));
        else
            DiffSign = sign(DiffDist);
        end
        % Same defect, same fix; see the splinefit branch above.
        pval = pkg.e_deviationpvalue(DiffDistAbs);

        T = table(gene, DiffDist, DiffDistAbs, DiffSign, pval);
        T = sortrows(T, 'DiffDistAbs', 'descend');

    otherwise
        error(['sc_dvg: unknown method ''%s''. Use ''splinefit'', ', ...
            '''analytic'' or ''brennecke''.'], method);
end
end

function [xyz, nearidx] = in_ordercurve(xyz, nearidx)
% Put foot points in curve order, so that plotting XYZ as a polyline traces
% the reference curve. The first coordinate is log1p(mu), strictly
% increasing in the curve parameter, so sorting on it is exactly curve
% order. NEARIDX is remapped so XYZ(NEARIDX,:) still returns each gene's
% own foot point, and the deviations scored from it are untouched.

[xyz, ord] = sortrows(xyz, 1);
backmap = zeros(numel(ord), 1);
backmap(ord) = 1:numel(ord);
nearidx = backmap(nearidx);
end

function [T, px1, py1, pz1, px2, py2, pz2] = ...
        in_scoredv(T1, T2, xyz1, xyz2, cL1, cL2, zeroatends, direction)
% Score two per-sample curve fits against each other and assemble the
% joined result table. T1 and T2 carry the columns SC_SPLINEFIT and
% SC_ANALYTICFIT both return, one row per gene, in the same gene order, and
% XYZ1/XYZ2 are the reference curves their NEARIDX indexes into.

px1 = T1.lgu; py1 = T1.lgcv; pz1 = T1.dropr;
px2 = T2.lgu; py2 = T2.lgcv; pz2 = T2.dropr;

v1 = ([px1 py1 pz1] - xyz1(T1.nearidx, :));
v2 = ([px2 py2 pz2] - xyz2(T2.nearidx, :));

DiffDist = vecnorm(v1 - v2, 2, 2);
% DIRECTION 'mean': +1 when the gene's mean is higher in sample 1 (tested)
% than in sample 2 (baseline), so an "up-regulated" DV gene means what it
% does in DE. LGU is log1p of the mean of library-size normalized counts,
% computed identically in both samples, and log1p is monotone.
% DIRECTION 'deviation': +1 when the gene sits further from its curve in
% sample 1 than in sample 2. LGU is sparse when the counts are, hence FULL.
switch direction
    case 'mean'
        DiffSign = full(sign(px1 - px2));
    case 'deviation'
        DiffSign = full(sign(vecnorm(v1, 2, 2) - vecnorm(v2, 2, 2)));
    otherwise
        % Unreachable: SC_DVG validates DIRECTION.
        error('sc_dvg:UnknownDirection', 'Unknown direction ''%s''.', direction);
end

% DIFFDIST is a norm, so it is non-negative and its null is the folded half
% of a symmetric distribution. This used to be zscore(DiffDist) read
% against a standard Normal, which is not a test at all: standardising a
% statistic by its own mean and SD makes the p-value a monotone function of
% the ranking, so the fraction called significant is fixed by the shape of
% DIFFDIST rather than by whether anything differs. On one homogeneous
% population split at random it called MORE genes than a real contrast did.
pval = pkg.e_deviationpvalue(DiffDist);

T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s', cL1{1}));
T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s', cL2{1}));
idxx = T1.(8) == 1 | T2.(8) == 1 | T1.(8) == max(T1.(8)) | T2.(8) == max(T2.(8));
T1 = T1(:, 1:end-1);
T2 = T2(:, 2:end);
T2 = T2(:, 1:end-1);
T1.Properties.VariableNames{1} = 'gene';

T = [T1 T2 table(DiffDist) table(DiffSign) table(pval)];

% idxx marks genes whose nearest point on either curve is its first or last
% point. For a spline that is the end of the fitted range, where the fit is
% not trustworthy, so the distance is discarded. DiffSign has to go with
% it: a gene with DiffDist 0 has no magnitude and therefore no direction,
% and leaving the sign behind let callers that split on sign alone report a
% zero-effect gene as up- or down-variable. That was 218 of 7,550 genes on
% the measured case, every one of them reported as a hit.
%
% The analytic curve has no such region -- it is defined at every mean --
% so ZEROATENDS is false there and every gene keeps its score.
if zeroatends
    T.DiffDist(idxx) = 0;
    T.DiffSign(idxx) = 0;
end
T = sortrows(T, 'DiffDist', 'descend');
end
