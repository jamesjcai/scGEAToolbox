function [T, info] = detect(X, genelist, opts)
%GLY.DETECT  Is this dataset deep enough for glyco analysis?
%
%   [T, info] = GLY.DETECT(X, genelist) reports how many glycogenes are
%   detected per cell, how that number scales with sequencing depth, and how
%   far the data sits from the point where deeper sequencing stops helping.
%   Run it before GLY.STATE or GLY:ENRICH: both rest on detection,
%   and neither can tell you whether the detection you have is the biology
%   or the library.
%
%   USAGE:
%     info = gly.detect(sce.X, sce.g);            % verdict for a dataset
%     [T, info] = gly.detect(sce.X, sce.g, Group=sce.c_cell_type_tx);
%     gly.detect(sce.X, sce.g, Plot=true);        % the depth curve
%
%   WHY IT MATTERS HERE MORE THAN ELSEWHERE. Glycogenes are lowly expressed,
%   so a glyco readout is further down the detection curve than a typical
%   marker-gene analysis and moves more for a given change in depth.
%   Measured on GSE218310: the two cultured lines sit at ~21000 median UMI
%   with ~110 glycogenes detected per cell, and the PBMC run at 4612 with 34.
%   Every module therefore reads as depleted in PBMC, and most of that gap is
%   the library rather than the glycobiology. That is the failure this
%   function exists to catch before it becomes a result.
%
%   THE SATURATION FIT IS DESCRIPTIVE, NOT MECHANISTIC. Binned median
%   detection is fit to D(u) = Dmax*u/(K+u), a saturating curve with no
%   claim about the sampling process behind it. Dmax estimates how many
%   glycogenes the cells could show at unlimited depth and K the depth at
%   which half of them are seen. What the fit is for is the ratio
%   info.saturation = D(median depth)/Dmax - the fraction of the reachable
%   repertoire the run actually captured. Read that, not the parameters.
%   Chrysinas et al. (NAR Genomics and Bioinformatics 2024, 6:lqae169)
%   report the same shape on Tabula Sapiens, with detection levelling off
%   near 220 glycogenes and little gained past ~65000 UMI per cell.
%
%   INPUTS:
%     X        - genes-by-cells counts. Only nonzero/zero is used for
%                detection; the row sums are used for depth, so pass RAW
%                counts. Library-size normalized input makes every cell the
%                same depth and the curve meaningless.
%     genelist - G-by-1 gene symbols
%     opts.Group     ([])   per-cell group label; adds a per-group table and
%                    is the quickest way to see a depth imbalance that would
%                    confound GLY.ENRICH
%     opts.GeneSets  ([])   table with Name and Genes columns, as
%                    GLY.GENESETS returns. Default is the curated
%                    glycobiology collection; the union of its genes is used.
%     opts.NumBins   (20)   depth bins for the curve
%     opts.MinBinCells (20) bins smaller than this are dropped from the fit
%     opts.Plot      (false) draw the curve
%
%   OUTPUT T: one row per depth bin - nCells, median UMI, median glycogenes
%     detected, median genes detected overall, and the glyco fraction.
%
%   OUTPUT info:
%     .nGlycoInData   glycogenes of the collection present in GENELIST
%     .medianUMI, .medianDetected   per cell, over all cells
%     .Dmax, .K       saturation fit parameters (NaN if the fit failed)
%     .saturation     fraction of Dmax reached at the median depth
%     .depthForHalf, .depthFor90   UMI needed for 50% and 90% of Dmax
%     .byGroup        per-group table when Group was given
%     .depthRatioUMI, .depthRatioDetected   max/min across groups. They can
%                     disagree sharply, and the DETECTED one is what governs
%                     whether a detection-based comparison is confounded:
%                     detection saturates, so groups six-fold apart in UMI
%                     can be a tenth of that apart in genes detected.
%     .verdict        one sentence, safe to print
%
% see also: GLY.ENRICH, GLY.STATE, GLY.GENESETS

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    opts.Group = []
    opts.GeneSets = []
    opts.NumBins (1, 1) double {mustBePositive} = 20
    opts.MinBinCells (1, 1) double {mustBePositive} = 20
    opts.Plot (1, 1) logical = false
end

if numel(genelist) ~= size(X, 1)
    error("GLY:DETECT:GeneCount", ...
        "GENELIST length (%d) must equal the number of rows of X (%d).", ...
        numel(genelist), size(X, 1));
end

glycogenes = i_glycogenes(opts.GeneSets);
isg = ismember(upper(genelist), upper(glycogenes));
nInData = sum(isg);
if nInData == 0
    error("GLY:DETECT:NoGlycogenes", ...
        "None of the %d glycogenes are in GENELIST.", numel(glycogenes));
end

% Full and double throughout: SCE.X is single sparse from R2025a on, and
% sparse medians and sparse logical masks are both more trouble than the
% memory they save on a vector this size.
umi = full(double(sum(X, 1)))';
nAll = full(double(sum(X > 0, 1)))';
nGly = full(double(sum(X(isg, :) > 0, 1)))';

% Relative, not absolute: library-size normalization leaves every column
% summing to the same 1e4-scale number up to float error, which is ~1e-11
% apart and would sail past an EPS comparison.
if (max(umi) - min(umi)) / max(max(umi), eps) < 1e-6
    warning("GLY:DETECT:UniformDepth", ...
        ['Every cell has the same library size, so X looks normalized. ', ...
        'Detection against depth is meaningless on normalized input; ', ...
        'pass raw counts.']);
end

% Equal-count depth bins, so every point on the curve rests on the same
% number of cells rather than on whatever happened to land in a fixed width.
edges = unique(quantile(umi, linspace(0, 1, opts.NumBins + 1)));
edges(1) = -Inf; edges(end) = Inf;
bin = discretize(umi, edges);
nb = max(bin);

binCells = zeros(nb, 1); binUMI = NaN(nb, 1);
binGly = NaN(nb, 1); binAll = NaN(nb, 1);
for k = 1:nb
    m = bin == k;
    binCells(k) = sum(m);
    if binCells(k) == 0, continue, end
    binUMI(k) = median(umi(m));
    binGly(k) = median(nGly(m));
    binAll(k) = median(nAll(m));
end
T = table((1:nb)', binCells, binUMI, binGly, binAll, binGly ./ max(binAll, 1), ...
    VariableNames = ["bin", "nCells", "medianUMI", "medianGlycoDetected", ...
    "medianGenesDetected", "glycoFraction"]);

use = binCells >= opts.MinBinCells & isfinite(binUMI) & isfinite(binGly);
if sum(use) < 3
    % Silently returning NaN here is the wrong failure: the caller asked a
    % question and would get no answer and no reason. NumBins defaults to 20
    % and MinBinCells to 20, so anything under ~400 cells lands here.
    warning("GLY:DETECT:TooFewBins", ...
        ['Only %d depth bin(s) hold at least %d cells, so the saturation ', ...
        'curve was not fit and Dmax is NaN. With %d cells, try ', ...
        'NumBins=%d or a lower MinBinCells.'], ...
        sum(use), opts.MinBinCells, numel(umi), ...
        max(3, floor(numel(umi) / opts.MinBinCells)));
end
[Dmax, K] = i_fitsaturation(binUMI(use), binGly(use), nInData);

medUMI = median(umi);
medDet = median(nGly);
saturation = NaN; d50 = NaN; d90 = NaN;
if isfinite(Dmax) && isfinite(K) && Dmax > 0
    saturation = (Dmax * medUMI / (K + medUMI)) / Dmax;
    d50 = K;                 % D(K) = Dmax/2 by construction
    d90 = 9 * K;             % D(9K) = 0.9*Dmax
end

info = struct();
info.nGlycoInData = nInData;
info.nGlycoInCollection = numel(glycogenes);
info.medianUMI = medUMI;
info.maxUMI = max(umi);
info.medianDetected = medDet;
info.Dmax = Dmax;
info.K = K;
info.saturation = saturation;
info.depthForHalf = d50;
info.depthFor90 = d90;
info.byGroup = table.empty;

if ~isempty(opts.Group)
    if numel(opts.Group) ~= size(X, 2)
        error("GLY:DETECT:GroupCount", ...
            "GROUP length (%d) must equal the number of columns of X (%d).", ...
            numel(opts.Group), size(X, 2));
    end
    [gi, gn] = findgroups(string(opts.Group(:)));
    info.byGroup = table(gn, accumarray(gi, 1), ...
        accumarray(gi, umi, [], @median), ...
        accumarray(gi, nGly, [], @median), ...
        accumarray(gi, nAll, [], @median), ...
        VariableNames = ["group", "nCells", "medianUMI", ...
        "medianGlycoDetected", "medianGenesDetected"]);
    % Two ratios, because they can disagree sharply and the detection one
    % is the one that matters. Detection saturates, so groups 6-fold apart
    % in UMI can be 1.1-fold apart in genes detected and barely confounded
    % at all; the reverse does not happen. GLY.ENRICH's tripwire is on
    % detection for the same reason.
    info.depthRatioUMI = max(info.byGroup.medianUMI) / ...
        max(min(info.byGroup.medianUMI), eps);
    info.depthRatioDetected = max(info.byGroup.medianGenesDetected) / ...
        max(min(info.byGroup.medianGenesDetected), eps);
end

info.verdict = i_verdict(info);

if opts.Plot
    i_plotcurve(umi, nGly, T(use, :), Dmax, K, info);
end

end


% ----------------------------------------------------------------------
function g = i_glycogenes(G)
if isempty(G)
    [setmatrx, ~, setgenes] = gly.genesets();
    g = setgenes(any(setmatrx, 1));
    return;
end
if ~istable(G) || ~all(ismember(["Name", "Genes"], G.Properties.VariableNames))
    error("GLY:DETECT:BadGeneSets", ...
        "GeneSets must be a table with Name and Genes columns.");
end
g = strings(0, 1);
for k = 1:height(G)
    v = strtrim(split(string(G.Genes(k)), ","));
    g = [g; v(strlength(v) > 0)]; %#ok<AGROW>
end
g = unique(g);
end


function [Dmax, K] = i_fitsaturation(u, d, nInData)
% D(u) = Dmax*u/(K+u), fit on binned medians by FMINSEARCH so that nothing
% outside base MATLAB is needed. Seeded from a Lineweaver-Burk style linear
% solve, which is biased but only ever used as a starting point.
Dmax = NaN; K = NaN;
if numel(u) < 3, return, end
u = u(:); d = d(:);
ok = u > 0 & d > 0;
u = u(ok); d = d(ok);
if numel(u) < 3, return, end

% A curve cannot be fit through points that all sit at one depth. This is
% what normalized input looks like by the time it reaches here, and without
% the guard the seeding solve below is rank deficient and warns.
if (max(u) - min(u)) / max(u) < 1e-6
    return;
end

A = [1 ./ u, ones(numel(u), 1)];
b = A \ (1 ./ d);
D0 = 1 / max(b(2), eps);
K0 = max(b(1) * D0, eps);
if ~isfinite(D0) || D0 <= 0, D0 = max(d); end
if ~isfinite(K0) || K0 <= 0, K0 = median(u); end

f = @(p) sum((d - abs(p(1)) .* u ./ (abs(p(2)) + u)).^2);
p = fminsearch(f, [D0, K0], optimset('Display', 'off', 'MaxIter', 2000));
Dmax = abs(p(1));
K = abs(p(2));

% A fit that runs away above the number of glycogenes actually present is
% extrapolating past anything the data could show, and its ratio would be
% meaningless. Report nothing rather than something unusable.
if ~isfinite(Dmax) || ~isfinite(K) || Dmax > 3 * nInData
    Dmax = NaN; K = NaN;
end
end


function s = i_verdict(info)
if ~isfinite(info.saturation)
    s = sprintf(['%d glycogenes present, median %d detected per cell at ', ...
        '%d UMI. The saturation curve could not be fit, so there is no ', ...
        'estimate of how much deeper sequencing would buy.'], ...
        info.nGlycoInData, round(info.medianDetected), round(info.medianUMI));
    return;
end
pct = round(100 * info.saturation);
% D(9K) is where 90% of the ceiling is reached, but on real data that lands
% an order of magnitude past the deepest cell observed, so it is the curve
% being extrapolated rather than a sequencing target anyone could act on.
% Say which it is instead of printing a number that reads like advice.
far = isfinite(info.depthFor90) && info.depthFor90 > 2 * info.maxUMI;
if pct >= 80
    tail = 'Deeper sequencing would add little; detection is near its ceiling.';
elseif pct >= 50 && ~far
    tail = sprintf(['Roughly %d UMI per cell would reach 90%% of the ', ...
        'reachable repertoire.'], round(info.depthFor90));
elseif pct >= 50
    tail = ['Saturation lies well beyond the depth of any cell here, so ', ...
        'treat detection-based comparisons between unequally sequenced ', ...
        'groups with care.'];
else
    tail = ['This is well short of saturation, so detection-based ', ...
        'comparisons are reporting depth as much as biology. Compare only ', ...
        'groups of matched depth, and check info.rhoDepth from ', ...
        'GLY.ENRICH before reading any enrichment as biology.'];
end
s = sprintf(['%d glycogenes present, median %d detected per cell at %d ', ...
    'UMI, which is %d%% of the estimated ceiling of %d. %s'], ...
    info.nGlycoInData, round(info.medianDetected), round(info.medianUMI), ...
    pct, round(info.Dmax), tail);
end


function i_plotcurve(umi, nGly, T, Dmax, K, info)
% Created here and nowhere earlier: a function that opens a window it never
% draws in is exactly the bug this toolbox just fixed in SC_KNNGRAPH.
f = figure('Name', 'gly.detect', 'Color', 'w', 'NumberTitle', 'off');
ax = axes(f);
n = numel(umi);
idx = 1:max(1, round(n / 5000)):n;          % thin for rendering only
scatter(ax, umi(idx), nGly(idx), 4, [0.7 0.75 0.8], 'filled', ...
    MarkerFaceAlpha = 0.35);
hold(ax, 'on');
plot(ax, T.medianUMI, T.medianGlycoDetected, 'o-', Color = [0.1 0.3 0.6], ...
    LineWidth = 1.5, MarkerFaceColor = 'w');
if isfinite(Dmax)
    uu = linspace(min(umi), max(umi), 200);
    plot(ax, uu, Dmax * uu ./ (K + uu), 'r-', LineWidth = 1.5);
    yline(ax, Dmax, 'r:', sprintf('ceiling %.0f', Dmax));
end
xline(ax, info.medianUMI, 'k--', 'median depth');
hold(ax, 'off');
set(ax, XScale = 'log');
xlabel(ax, 'UMI per cell');
ylabel(ax, 'glycogenes detected');
title(ax, 'Glycogene detection against sequencing depth');
subtitle(ax, info.verdict, FontWeight = 'normal', FontSize = 8);
grid(ax, 'on');
end
