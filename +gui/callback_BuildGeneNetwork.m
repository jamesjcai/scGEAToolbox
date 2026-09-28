function callback_BuildGeneNetwork(src, ~)
% CALLBACK_BUILDGENENETWORK  Build a GRN for one cell population.
%
% One flow for Network > Build Gene Regulatory Network (GRN)...:
%
%   1. Cells   all cells, or the cells in chosen groups of one grouping
%   2. Genes   pick from a list, paste names, or all genes
%   3. Method  selected genes: any method in net.grnmethods, then the
%              transform; gui.i_showgrn saves the network to the
%              working folder and draws it.
%              All genes: gui.callback_BuildGRNAllGenes takes over.
%
% See also gui.callback_BuildGRNAllGenes, gui.callback_CompareGeneNetwork,
%   net.grnmethods, gui.i_showgrn, sc_grnview.

[FigureHandle, sce] = gui.gui_getfigsce(src);

% --- Cells
cellIdx = gui.i_pickgrncells(sce, FigureHandle, sprintf(['A gene ', ...
    'regulatory network describes one cell population. Build it from ', ...
    'all %d cells, or only from cells in groups you select?'], sce.NumCells));
if isempty(cellIdx), return; end

% --- Genes
[glist, useAll] = gui.i_pickgrngenes(sce.g, FigureHandle);
if useAll
    gui.callback_BuildGRNAllGenes(src, cellIdx);
    return;
end
if isempty(glist), return; end

[y, i] = ismember(upper(glist), upper(sce.g));
if ~all(y), error('Runtime error.'); end
fprintf("%s\n", glist)

% The method list, its preferred transform and its pre-run warnings all
% live in net.grnmethods.
method = gui.i_selectgrnmethod(FigureHandle, sprintf(['How should the ', ...
    'links between the %d genes be scored?'], numel(glist)));
if isempty(method), return; end

[Xt] = gui.i_transformx(sce.X(:, cellIdx), true, method.Transform, FigureHandle);
if isempty(Xt), return; end
x = Xt(i, :);
if ~gui.i_confirmgrnrun(method, x, 1, FigureHandle), return; end

fw = gui.myWaitbar(FigureHandle);
try
    A = sc_grn(x, method.Key);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

% Saved silently to the working folder, then drawn whole (or sketched,
% for a long gene list); the figure's toolbar exports and locates it.
key = char(method.Key);
try
    gui.i_showgrn(A, glist, key, lower(matlab.lang.makeValidName(key)), ...
        FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end
end
