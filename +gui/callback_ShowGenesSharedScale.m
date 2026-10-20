function callback_ShowGenesSharedScale(src, ~)
%CALLBACK_SHOWGENESSHAREDSCALE Show genes side by side on one colour scale.
%
%   Behind View > Show Gene Expression > Genes on One Colour Scale. Two to
%   four genes are drawn on the embedding side by side, all on the same
%   colour range with one colour bar, so one colour means one expression
%   level in every panel. Selected Genes scales each gene to its own range
%   instead, which makes a weakly expressed gene look as strong as a highly
%   expressed one.
%
%   Zooming, panning or rotating one panel moves them all. Toolbar buttons
%   switch to a range per gene and back, and set the shared range by hand.
%   With two genes the title says how many cells detect each, and both.
%
%   The idea is AnnZarro's Lock Range (Otto, Baasri and Setty, 2026), which
%   holds two panels on one fixed colour range to compare genes.
%
%   See also GUI.I_SHAREDSCALEPLOT, GUI.CALLBACK_SHOWGENEEXPR.

[FigureHandle, sce] = gui.gui_getfigsce(src);
minGenes = 2;
maxGenes = 4;
if isempty(gui.i_checkexistingembed(sce))
    gui.myWarndlg(FigureHandle, ['The cells have no embedding yet. ' ...
        'Embed them first, then try again.']);
    return;
end

viewAngle = [];
mainAxes = findall(FigureHandle, 'type', 'axes');
if ~isempty(mainAxes)
    [az, el] = view(mainAxes(1));
    viewAngle = [az, el];
end

glist = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    return;
end
glist = unique(string(glist), 'stable');
if numel(glist) < minGenes || numel(glist) > maxGenes
    gui.myWarndlg(FigureHandle, sprintf(['Pick %d to %d genes to compare; ' ...
        '%d were picked.'], minGenes, maxGenes, numel(glist)));
    return;
end

[Xt, used] = gui.i_transformx(sce.X, [], [], FigureHandle);
if isempty(Xt)
    return;
end
numGenes = numel(glist);
values = zeros(sce.NumCells, numGenes);
for k = 1:numGenes
    values(:, k) = full(Xt(find(sce.g == glist(k), 1), :)).';
end
if used == "none"
    valueLabel = 'Expression (data unchanged)';
else
    valueLabel = 'Expression (transformed)';
end

gui.i_sharedscaleplot(values, glist, sce.s, FigureHandle, viewAngle, valueLabel);
end
