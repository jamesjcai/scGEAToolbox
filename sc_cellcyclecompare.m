function [result, fig] = sc_cellcyclecompare(Xcontrol, gControl, Xtreated, gTreated, options)
%SC_CELLCYCLECOMPARE  Compare cell cycle position between two samples.
%   result = SC_CELLCYCLECOMPARE(Xcontrol, gControl, Xtreated, gTreated)
%   places every cell of both samples on the cell cycle with tricycle
%   (run.ml_tricycle), using the control as the reference, and tests whether
%   the treated cells sit somewhere else on the cycle.
%
%   [result, fig] = SC_CELLCYCLECOMPARE(...) also returns the figure.
%
%   Xcontrol, Xtreated : raw count matrices, genes x cells.
%   gControl, gTreated : their gene names. Only genes in both are used,
%                        matched exactly and then ignoring case.
%
%   Name-value arguments:
%   Species         : "human" or "mouse"; guessed from the symbols if omitted.
%   NumPermutations : permutations for the p-values (default 2000).
%   GroupNames      : names for the two samples (default ["Control", "Treated"]).
%   Plot            : draw the figure (default true).
%   Parent          : figure to draw into, e.g. a gui.myFigure window's
%                     FigHandle (default: a new figure).
%
%   Both samples go through tricycle in one run: normalised together, and
%   with the angle taken around the control cells' mean. tricycle on its own
%   centres on whatever cells it is given, so scoring the treated sample by
%   itself would move its origin with the treatment.
%
%   RESULT fields, each two-element row in [control, treated] order:
%   GroupNames, NumCells, NumPermutations
%   Position            : 1x2 cell of positions in radians, [0, 2*pi)
%   CircularMean        : mean direction of the positions
%   Concentration       : mean resultant length, 0 (spread round the cycle)
%                         to 1 (all at one position)
%   MeanShift           : angle between the two mean directions, [0, pi]
%   PValueMeanShift     : permutation p-value for MeanShift
%   PValueConcentration : permutation p-value for the difference in
%                         Concentration -- a partial arrest can bunch cells
%                         up without moving their mean
%   BinEdges, BinShares : share of cells in eight bins of pi/4 (2 x 8)
%
%   Positions are compared, not phase labels: the label boundaries are
%   approximate (see RUN.ML_TRICYCLE). The figure marks where the S and G2M
%   labels fall only as a guide.
%
%   See also RUN.ML_TRICYCLE, SC_CELLCYCLESCORE

arguments
    Xcontrol
    gControl
    Xtreated
    gTreated
    options.Species {mustBeTextScalar} = ""
    options.NumPermutations (1, 1) {mustBeInteger, mustBePositive} = 2000
    options.GroupNames (1, 2) string = ["Control", "Treated"]
    options.Plot (1, 1) logical = true
    options.Parent = []
end

[Xcontrol, Xtreated, genelist] = i_commongenes(Xcontrol, gControl, Xtreated, gTreated);
nControl = size(Xcontrol, 2);

% The control is in the sample and also is the reference, so both groups
% come back on one normalisation and one origin.
[~, position] = run.ml_tricycle([Xcontrol, Xtreated], genelist, ...
    Species=options.Species, ReferenceX=Xcontrol, ReferenceGenes=genelist);
posControl = position(1:nControl);
posTreated = position(nControl+1:end);

result = i_compare(posControl, posTreated, options.NumPermutations);
result.GroupNames = options.GroupNames;

fig = [];
if options.Plot
    fig = i_plot(result, options.Parent);
end
end

function result = i_compare(posA, posB, numPermutations)
zA = exp(1i*posA(:));
zB = exp(1i*posB(:));
nA = numel(zA);
[shift, concDiff] = i_stats(zA, zB);

% Permutation null: shuffle which sample each cell came from.
zAll = [zA; zB];
nullShift = zeros(numPermutations, 1);
nullConc = zeros(numPermutations, 1);
for k = 1:numPermutations
    q = zAll(randperm(numel(zAll)));
    [nullShift(k), nullConc(k)] = i_stats(q(1:nA), q(nA+1:end));
end

edges = (0:8)*pi/4;
result = struct( ...
    "GroupNames", strings(1, 2), ...
    "NumPermutations", numPermutations, ...
    "NumCells", [numel(zA), numel(zB)], ...
    "Position", {{posA(:), posB(:)}}, ...
    "CircularMean", mod(angle([mean(zA), mean(zB)]), 2*pi), ...
    "Concentration", abs([mean(zA), mean(zB)]), ...
    "MeanShift", shift, ...
    "PValueMeanShift", (1 + sum(nullShift >= shift))/(numPermutations + 1), ...
    "PValueConcentration", (1 + sum(nullConc >= concDiff))/(numPermutations + 1), ...
    "BinEdges", edges, ...
    "BinShares", [histcounts(posA, edges, Normalization="probability"); ...
                  histcounts(posB, edges, Normalization="probability")]);
end

function [shift, concDiff] = i_stats(zA, zB)
mA = mean(zA);
mB = mean(zB);
shift = abs(angle(mA*conj(mB)));
concDiff = abs(abs(mA) - abs(mB));
end

function fig = i_plot(result, fig)
if isempty(fig)
    fig = figure;
end
fig.Name = "Cell cycle position";
% White, as the toolbox's other figures read: the two-panel layout leaves
% much of the figure's default light gray showing.
fig.Color = "w";
layout = tiledlayout(fig, 1, 2, TileSpacing="compact");
% Two lines: one does not fit a window of the default width.
title(layout, ["Cell cycle position (tricycle)", sprintf("mean shift %.2f\\pi (%s); " + ...
    "concentration %.2f vs %.2f (%s)", result.MeanShift/pi, ...
    i_formatp(result.PValueMeanShift, result.NumPermutations), ...
    result.Concentration(1), result.Concentration(2), ...
    i_formatp(result.PValueConcentration, result.NumPermutations))]);
names = result.GroupNames;
plotEdges = (0:24)*pi/12;
colors = lines(2);

pax = polaraxes(layout);
pax.Layout.Tile = 1;
hold(pax, "on");
for k = 1:2
    polarhistogram(pax, result.Position{k}, plotEdges, Normalization="probability", ...
        FaceColor=colors(k, :), FaceAlpha=0.45, EdgeColor="none");
end
rMax = max(pax.RLim);
for k = 1:2
    % Mean direction, drawn to the concentration's share of the radius.
    polarplot(pax, [0, result.CircularMean(k)], [0, result.Concentration(k)*rMax], ...
        Color=colors(k, :), LineWidth=2);
end
pax.ThetaAxisUnits = "radians";
legend(pax, names, Location="southoutside");

ax = nexttile(layout, 2);
hold(ax, "on");
for k = 1:2
    histogram(ax, result.Position{k}, plotEdges, Normalization="probability", ...
        DisplayStyle="stairs", EdgeColor=colors(k, :), LineWidth=1.5, ...
        DisplayName=sprintf("%s (n = %d)", names(k), result.NumCells(k)));
end
xlim(ax, [0, 2*pi]);
i_shade(ax, [0.5, 1]*pi, "S", 0.12);
i_shade(ax, [1, 1.25]*pi, "G2M", 0.3);
xticks(ax, (0:4)*pi/2);
xticklabels(ax, ["0", "\pi/2", "\pi", "3\pi/2", "2\pi"]);
xlabel(ax, "Position (radians)");
ylabel(ax, "Share of cells");
legend(ax, Location="northoutside");
end

function txt = i_formatp(p, numPermutations)
% At the floor no permutation reached the observed value, so the p-value is
% bounded by the number of permutations rather than measured.
if p <= 1/(numPermutations + 1)
    txt = sprintf("p \\leq %.2g", p);
else
    txt = sprintf("p = %.2g", p);
end
end

function i_shade(ax, range, label, alpha)
% Where run.ml_tricycle's labels fall; a guide only.
xregion(ax, range(1), range(2), FaceColor=[0.5, 0.5, 0.5], FaceAlpha=alpha, ...
    HandleVisibility="off");
yl = ylim(ax);
text(ax, mean(range), yl(2), label, HorizontalAlignment="center", ...
    VerticalAlignment="top", Color=[0.35, 0.35, 0.35]);
end

function [Xa, Xb, genelist] = i_commongenes(Xa, ga, Xb, gb)
ga = string(ga(:));
gb = string(gb(:));
if size(Xa, 1) ~= numel(ga) || size(Xb, 1) ~= numel(gb)
    error("sc_cellcyclecompare:size", ...
        "Each count matrix needs one gene name per row: got %d/%d (control) and %d/%d (treated).", ...
        size(Xa, 1), numel(ga), size(Xb, 1), numel(gb));
end
% Exact names first, so two genes that differ only in case stay apart;
% then ignoring case.
[found, loc] = ismember(ga, gb);
[foundCI, locCI] = ismember(upper(ga(~found)), upper(gb));
loc(~found) = locCI;
found(~found) = foundCI;
if ~any(found)
    error("sc_cellcyclecompare:noCommonGenes", ...
        "The two samples share no gene names. Check that both use the same gene symbols.");
end
Xa = Xa(found, :);
Xb = Xb(loc(found), :);
genelist = ga(found);
end
