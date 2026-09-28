function [T, X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2] = ...
        sc_dvg(sce1, sce2, cL1, cL2, method, direction, options)
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
% Name-value options:
%   NumPermutations : 0 (default) or at least 10; 'splinefit' and
%                'analytic' only. When positive, the pval column comes from
%                a permutation null instead of the closed-form one;
%                DiffDist, DiffSign and the ranking are unchanged. See
%                "Calibration" below.
%
% Outputs:
%   T    : results table sorted by DiffDist descending. DiffDist is the
%          size of the variability difference; DiffSign is its direction
%          as chosen by DIRECTION. Split DV genes on DiffSign.
%   X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2 :
%          curve visualization data ('splinefit' and 'analytic'; empty for
%          'brennecke')
%
% Calibration. With the default NumPermutations=0, pval reads every gene's
% distance against one half-normal null whose scale is the median over all
% genes (PKG.E_DEVIATIONPVALUE). Real genes do not share one noise scale:
% it rises with expression, and genes carried by a few contaminating cells
% are heavy-tailed. On random halves of one cell type in the bundled
% example data, where nothing differs, 5.8-12.4% of genes came back at
% p < 0.05 and 63-383 passed BH q < 0.05, across all three methods. Treat
% those p-values as a ranking aid, not a test.
%
% NumPermutations=K reassigns the pooled cells to two groups of the
% original sizes K times, reruns the same method, and records each gene's
% distance. Each gene is scaled by the RMS of its own permuted distances,
% and its scaled distance is read against the scaled permuted distances of
% all genes pooled. On the same random halves, K=30 gave 4.6-5.5% at
% p < 0.05 and no BH calls for 'splinefit' and 'analytic'. The cost is
% K+1 runs of the method, about 30x at K=30. The permutations draw from
% their own seeded stream, so results are reproducible and the global
% random state is not touched.
%
% What it tests is the statistic as defined: whether a gene's deviation
% from its own group's curve differs more than relabelling cells explains.
% When many genes change, the curve itself moves, and an unchanged gene's
% deviation from it moves too -- measurably so for well-expressed genes,
% whose own noise is small. The permutation p-value calls those; the
% default's single wide scale absorbs them. On Poisson fixtures with
% variance planted in 5-20% of genes (fixtures on which the default is
% calibrated), BH precision was 0.79-0.98 with permutations against
% 0.91-1.00 without, and recall 0.90-0.97 against 0.96-1.00.
%
% 'brennecke' is excluded: its residual CV^2 is taken against a trend
% fitted to all genes together, so an effect in some genes moves every
% gene's residual, and pooling the cells removes that shift. On the
% planted effect the permutation p-value called 1688 of 2000 genes at
% precision 0.22.

arguments
    sce1
    sce2
    cL1 = []
    cL2 = []
    method = []
    direction = []
    options.NumPermutations (1,1) double {mustBeInteger, mustBeNonnegative} = 0
end

if nargin < 3 || isempty(cL1), cL1 = {'1'}; end
if nargin < 4 || isempty(cL2), cL2 = {'2'}; end
if nargin < 5 || isempty(method), method = 'splinefit'; end
if nargin < 6 || isempty(direction), direction = 'mean'; end
direction = validatestring(direction, {'mean', 'deviation'}, ...
    mfilename, 'direction', 6);
minPermutations = 10;
if options.NumPermutations > 0 && options.NumPermutations < minPermutations
    error('sc_dvg:TooFewPermutations', ...
        ['NumPermutations must be 0 or at least %d; each gene''s null scale ', ...
        'is estimated from its own permuted distances.'], minPermutations);
end
if options.NumPermutations > 0 && strcmpi(method, 'brennecke')
    error('sc_dvg:PermutationsNotSupported', ...
        ['NumPermutations is not supported for ''brennecke''. Use ''splinefit'' ', ...
        'or ''analytic'', or set NumPermutations to 0.']);
end

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

[T, X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2] = ...
    in_dvgcore(X1_ori, X2_ori, g_ori, cL1, cL2, method, direction);
if options.NumPermutations > 0
    T.pval = in_permutationpval(T, X1_ori, X2_ori, g_ori, method, direction, ...
        options.NumPermutations);
end
end

function [T, X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2] = ...
        in_dvgcore(X1_ori, X2_ori, g_ori, cL1, cL2, method, direction)
% Score one pair of raw count matrices (genes x cells, rows matching
% G_ORI) with METHOD. Everything SC_DVG returns comes from here, and the
% permutation null reruns exactly this on relabelled cells.

X1 = []; X2 = []; g = []; xyz1 = []; xyz2 = [];
px1 = []; py1 = []; pz1 = []; px2 = []; py2 = []; pz2 = [];

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

        % SC_SPLINEFIT drops the genes that are all-zero in its own sample,
        % so a gene detected in only one group leaves G1 and G2 different
        % lengths. Keep the genes both fits scored; NEARIDX indexes each
        % fit's own XYZ, so subsetting rows leaves it valid.
        [~, ia, ib] = intersect(g1, g2, 'stable');
        T1 = T1(ia, :);
        X1 = X1(ia, :);
        T2 = T2(ib, :);
        X2 = X2(ib, :);
        g = g1(ia);

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

function pval = in_permutationpval(T, X1, X2, g, method, direction, numPerm)
% Permutation-studentized p-value for the distance column of T. See
% "Calibration" in the help. X1 and X2 are the raw counts the observed run
% used; their cells are pooled and reassigned to groups of the original
% sizes NUMPERM times.
%
% A distance of 0 is a gene the method discarded (the ends of the spline),
% so it is left out of that permutation's null and, when observed, gets
% p = 1. So does a gene with too few permuted distances to scale it.

if strcmpi(method, 'brennecke')
    magnitude = 'DiffDistAbs';
else
    magnitude = 'DiffDist';
end
observed = T.(magnitude);
Xall = [X1, X2];
numCells1 = size(X1, 2);
numCells = size(Xall, 2);
stream = RandStream('mt19937ar', 'Seed', 0);

D = nan(height(T), numPerm);
for r = 1:numPerm
    order = randperm(stream, numCells);
    Tperm = in_dvgcore(Xall(:, order(1:numCells1)), Xall(:, order(numCells1+1:end)), ...
        g, {'1'}, {'2'}, method, direction);
    [found, loc] = ismember(T.gene, Tperm.gene);
    d = nan(height(T), 1);
    d(found) = Tperm.(magnitude)(loc(found));
    d(d == 0) = NaN;
    D(:, r) = d;
end

% Each gene's scale is the RMS of its own permuted distances. The null
% values are scaled leaving their own permutation out, so that no value is
% divided by a scale it helped set.
sumSq = sum(D.^2, 2, 'omitnan');
count = sum(isfinite(D), 2);
minCount = max(5, ceil(numPerm/2));
usable = count >= minCount & sumSq > 0;
scale = sqrt(sumSq./count);
looScale = sqrt((sumSq - D.^2)./(count - 1));
nullZ = D(usable, :)./looScale(usable, :);
nullZ = nullZ(isfinite(nullZ));

pval = ones(height(T), 1);
scored = usable & observed > 0;
if isempty(nullZ) || ~any(scored)
    return;
end
z = observed(scored)./scale(scored);
% Share of pooled null values at or above each observed value, with the
% usual +1 so that no p-value is 0. The null values go first so that the
% stable sort counts a tie as at-or-above.
numNull = numel(nullZ);
[~, ord] = sort([nullZ; z], 'descend');
isNull = [true(numNull, 1); false(numel(z), 1)];
nullAbove = cumsum(isNull(ord));
atOrAbove = zeros(numel(z), 1);
atOrAbove(ord(~isNull(ord)) - numNull) = nullAbove(~isNull(ord));
pval(scored) = (1 + atOrAbove)/(numNull + 1);
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
%
% IDXX marks the genes whose distance is discarded below (see ZEROATENDS).
% They get p = 1, and the null is fitted without them. The p-value used to
% be computed before they were zeroed, so a gene reported with DiffDist 0
% and DiffSign 0 could carry one of the smallest p-values in the table and
% pass an FDR cut on PVAL; and the untrusted distances it discards were
% setting the null scale for every other gene.
idxx = T1.(8) == 1 | T2.(8) == 1 | T1.(8) == max(T1.(8)) | T2.(8) == max(T2.(8));
if zeroatends
    scored = ~idxx;
else
    scored = true(size(DiffDist));
end
pval = ones(size(DiffDist));
pval(scored) = pkg.e_deviationpvalue(DiffDist(scored));

T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s', cL1{1}));
T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s', cL2{1}));
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
