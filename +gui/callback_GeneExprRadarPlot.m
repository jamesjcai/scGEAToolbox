function callback_GeneExprRadarPlot(src, ~)
%CALLBACK_GENEEXPRRADARPLOT Radar plot of mean gene expression per cell group.
%
%   Behind Plots > Gene Expression Radar Plot. One axis per gene, one
%   polygon per cell group, each vertex the group's mean log-normalised
%   expression of that gene. The groups are picked with the same chooser
%   as the Dotplot, Heatmap and violin plots (gui.i_selectgroupsubset), and
%   the drawing, saving and renaming are gui.i_spiderplot's, shared with
%   the Gene Program Radar Plot.
%
%   See also GUI.I_SPIDERPLOT, GUI.CALLBACK_GETCELLSIGNATUREMATRIX.

minGenes = 3;   % fewer axes make no radar

[FigureHandle, sce] = gui.gui_getfigsce(src);

% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL").
[thisc, clabel] = gui.i_selectnclass(sce, false, [], [], FigureHandle);
if isempty(thisc), return; end
thisc = string(thisc);

[picked, levels] = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
if isempty(picked), return; end
thisc = thisc(picked);

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist), return; end
[~, gidx] = ismember(upper(glist), upper(sce.g));
glist = glist(gidx > 0);
gidx = gidx(gidx > 0);
if numel(gidx) < minGenes
    gui.myWarndlg(FigureHandle, sprintf(['A radar plot needs at least ' ...
        '%d genes, one per axis. Select %d or more genes.'], minGenes, minGenes));
    return;
end

[Xt] = gui.i_transformx(sce.X, true, "libsize_log1p", FigureHandle);
if isempty(Xt), return; end
% Subset after normalising, so each cell's values match an all-cells plot.
Y = full(Xt(gidx, picked)).';

% I_SPIDERPLOT reads only the cell IDs from its SCE argument, for the
% per-cell table it saves; a struct carrying them spares copying the
% count matrix just to subset it.
cellinfo = struct('c_cell_id', sce.c_cell_id(picked));

try
    gui.i_spiderplot(Y, thisc, cellstr(glist(:).'), cellinfo, ...
        FigureHandle, levels, ValueKind="expression");
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end
end
