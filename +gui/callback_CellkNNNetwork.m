function callback_CellkNNNetwork(src, ~)
%CALLBACK_CELLKNNNETWORK Draw the cells' k-nearest-neighbour graph over the embedding.
%
%   Behind Plots > Cell kNN Network. Each cell is joined to its k nearest
%   neighbours in the current embedding (every column of SCE.S), and the
%   graph is drawn at the cells' embedding positions.
%
%   Narrowing to a subset of cells is an opt-in first step, through the
%   group chooser the Dotplot, Heatmap and violin plots use (see
%   gui.i_selectgroupsubset). The graph is then built on the picked cells
%   alone, so each one's neighbours are its nearest picked cells.
%
%   See also SC_KNNGRAPH, GUI.I_SINGLEGRAPH.

[FigureHandle, sce] = gui.gui_getfigsce(src);

picked = true(sce.NumCells, 1);
answer = gui.myQuestdlg(FigureHandle, ...
    "Build the network on all cells, or only cells in selected groups?", "", ...
    {'All Cells', 'Selected Groups', 'Cancel'}, 'All Cells');
switch answer
    case 'All Cells'
        % Keep every cell.
    case 'Selected Groups'
        [thisc, clabel] = gui.i_selectnclass(sce, false, [], [], FigureHandle);
        if isempty(thisc), return; end
        picked = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
        if isempty(picked), return; end
    otherwise
        return;
end

k = gui.i_inputnumk(4, 0, [], [], FigureHandle);
if isempty(k), return; end
s = sce.s(picked, :);
if k >= size(s, 1)
    gui.myWarndlg(FigureHandle, sprintf(['k must be smaller than the ' ...
        'number of cells (%d). Pick more groups or a smaller k.'], size(s, 1)));
    return;
end

fw = gui.myWaitbar(FigureHandle);
try
    [A] = sc_knngraph(s, k, false, 1, FigureHandle);
    a = string(sce.c_cell_id(picked));
    a = strrep(a, '_', '\_');
    G = graph(A, a);
    p = gui.i_singlegraph(G, '', FigureHandle);
    p.XData = s(:, 1)';
    p.YData = s(:, 2)';
    % The graph is built on every column of S, so draw a 3-D embedding in
    % 3-D rather than flattening it.
    if size(s, 2) >= 3
        p.ZData = s(:, 3)';
        view(p.Parent, 3);
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
end
