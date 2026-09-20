function [expanded] = callback_ExpandBrushedCells(src)
%CALLBACK_EXPANDBRUSHEDCELLS Grow the brush to the whole group it landed in.
%
%   expanded = gui.callback_ExpandBrushedCells(src)
%
%   Brushing a handful of cells of a cluster and expanding the selection to
%   the rest of that cluster is what GUI.I_EXPANDBRUSHED offers inside the
%   handlers that consume a selection - Delete Selected Cells asks it every
%   time. This does the same thing on its own and writes the result back to
%   the plot's BrushData, so the expanded selection is what every later
%   brush-driven command sees, whichever one that turns out to be.
%
%   The groups are the current grouping, SCE.C, which is also what colours
%   the plot. A brush that spans several groups expands to all of them, after
%   confirming: the alternative, refusing, leaves the user to redo a selection
%   whose intent was not in doubt.
%
%   EXPANDED is the new selection as a logical column, or [] when nothing was
%   changed.
%
%   src can be the scgeatool app or a component of a plain figure.
%
%   See also gui.i_expandbrushed, gui.callback_SelectCellsByClass.

expanded = [];
[FigureHandle, sce] = gui.gui_getfigsce(src);

h = in_scatterhandle(src, FigureHandle);
if ~pkg.i_isvalid(h)
    gui.myWarndlg(FigureHandle, 'No plot available.');
    return;
end

brushed = logical(h.BrushData(:));
if ~any(brushed)
    gui.myHelpdlg(FigureHandle, ['No cells are highlighted. Use the data ' ...
        'brush tool to highlight a few cells of the group you want, then ' ...
        'run this again.'], 'Expand Highlighted Cells');
    return;
end

[c, cL] = findgroups(string(sce.c));
c = c(:);
if numel(c) ~= numel(brushed)
    % A plot of something other than the cells - a gene scatter, say. There
    % is no group to expand to, and indexing on would be silently wrong.
    gui.myWarndlg(FigureHandle, ['The highlighted points are not cells of ' ...
        'this dataset, so there is no cell group to expand to.'], ...
        'Expand Highlighted Cells');
    return;
end
if isscalar(unique(c))
    gui.myWarndlg(FigureHandle, ['Every cell is in the same group, so ' ...
        'expanding the selection would select everything. Cluster the ' ...
        'cells, or switch to another grouping variable, first.'], ...
        'Expand Highlighted Cells');
    return;
end

groups = unique(c(brushed));
mask = ismember(c, groups);

if isequal(mask, brushed)
    gui.myHelpdlg(FigureHandle, sprintf(['The highlighted cells are ' ...
        'already the whole of %s.'], in_namelist(cL, groups)), ...
        'Expand Highlighted Cells');
    return;
end

if ~isscalar(groups)
    % Say which groups before growing the selection by a few thousand cells.
    answer = gui.myQuestdlg(FigureHandle, sprintf(['The %d highlighted ' ...
        'cells fall in %d groups: %s. Expand the selection to all %d cells ' ...
        'of those groups?'], sum(brushed), numel(groups), ...
        in_namelist(cL, groups), sum(mask)), 'Expand Highlighted Cells', ...
        {'Expand', 'Cancel'}, 'Expand');
    if ~strcmp(answer, 'Expand'), return; end
end

% The brush is the selection every other handler reads, so writing it back
% here is what makes this worth having as a menu item of its own.
gui.i_setbrushdata(h, mask);
expanded = mask;

gui.myHelpdlg(FigureHandle, sprintf(['Selection expanded from %d to %d ' ...
    'cells - all of %s.'], sum(brushed), sum(mask), ...
    in_namelist(cL, groups)), 'Expand Highlighted Cells');
end

function [h] = in_scatterhandle(src, FigureHandle)
% The app keeps its scatter in APP.H; a plain figure has to be searched.

if isa(src, 'matlab.apps.AppBase') && isprop(src, 'h') && ...
        pkg.i_isvalid(src.h)
    h = src.h;
    return;
end
h = findobj(FigureHandle, 'Type', 'Scatter');
if ~isscalar(h) && ~isempty(h)
    h = h(1);
end
end

function [s] = in_namelist(cL, groups)
s = pkg.i_namesummary(cL(groups));
end
