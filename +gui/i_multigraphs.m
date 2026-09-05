function [hFig, pv] = i_multigraphs(Gs, nodeinfo, figname, parentfig)
% I_MULTIGRAPHS  Multi-panel network viewer sharing one node layout.
%
%   [hFig, pv] = gui.i_multigraphs(Gs, nodeinfo, figname, parentfig)
%
%   Gs        Cell array of graph/digraph objects that share an identical
%             Nodes.Name column. Usually {G1, G2} - two conditions - and
%             optionally a third holding their difference.
%   nodeinfo  Struct of per-node decoration, all fields optional:
%               .PanelTitles 1xP string, title above each axes
%               .Highlight   Nx1 logical, the nodes to pick out: red fill
%                            and a bold red label. Everything else is grey.
%               .HighlightName  legend text for the highlighted nodes
%               .PlainName   legend text for the rest
%               .IsTF        Nx1 logical, TFs get diamond markers
%               .EdgeRef     1xP numeric, the |weight| that maps to the
%                            widest line in each panel. Panels sharing a
%                            value are directly comparable.
%               .Legend      scalar string appended to the info dialog
%   figname   Figure title.
%   parentfig Parent figure to inherit the theme from and centre on.
%
%   Returns the figure handle and the 1xP array of GraphPlot handles.
%
% See also ten.sctenifoldnetview, gui.i_doublegraphs, gui.i_singlegraph.

if nargin < 4, parentfig = []; end
if nargin < 3, figname = ''; end
if nargin < 2, nodeinfo = struct; end
if nargin < 1, error('USAGE: gui.i_multigraphs(Gs, nodeinfo)'); end

if ~iscell(Gs), Gs = {Gs}; end
npanel = numel(Gs);
nodename = string(Gs{1}.Nodes.Name);
nnode = numel(nodename);
for k = 2:npanel
    if ~isequal(string(Gs{k}.Nodes.Name), nodename)
        error('All panels must share an identical node list.');
    end
end

nodeinfo = i_filldefaults(nodeinfo, npanel, nnode);

% Node style is fixed for the life of the figure; edge style is recomputed
% whenever a panel is redrawn, so it is derived on the fly below.
% One size for every node. Colour already says which nodes matter, and a
% second channel saying the same thing only adds a size difference the
% reader has to decide whether to read anything into.
markerSize = repmat(6, nnode, 1);
nodeMarker = repmat({'o'}, nnode, 1);
nodeMarker(nodeinfo.IsTF) = {'d'};

hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;
% After the figure exists: the accent red depends on its theme.
nodeColor = i_nodecolor(nodeinfo.Highlight, gui.i_accentred(hFig));
set(0, 'CurrentFigure', hFig);
hFig.Position(3) = hFig.Position(3) * min(npanel, 3) * 0.8;

t = tiledlayout(hFig, 1, npanel, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(1, npanel);
pv = gobjects(1, npanel);
for k = 1:npanel
    ax(k) = nexttile(t);
end

% One force layout on the union of all panels, so a node sits at the same
% spot everywhere and the eye can compare panels position by position.
xy = i_unionlayout(Gs, nodename, '');

labelMode = 1; % 1 = all, 2 = highlighted only, 3 = none
edgeWidthMax = 2;
selectedNode = 0;

for k = 1:npanel
    in_drawpanel(k);
end
title(t, figname, 'Interpreter', 'none');

lgd = in_buildlegend();

hx.addCustomButton('off', @in_ChangeFontSize, 'noun_font_size_591141.gif', 'Change Font Size of Nodes');
hx.addCustomButton('off', @in_ChangeWeight, 'weight_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Width of Edges');
hx.addCustomButton('off', @in_ChangeLabels, 'label.jpg', 'Label All Nodes / Highlighted Only / None');
hx.addCustomButton('off', @in_ChangeLayout, 'group_work_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Network Layout');
hx.addCustomButton('on', @in_ClearSelection, 'refresh.jpg', 'Clear Node Selection');
hx.addCustomButton('off', @in_SaveAdj, 'floppy-disk-arrow-in.jpg', 'Export & Save Data');
hx.addCustomButton('on', @in_ToggleLegend, 'icon-fa-list-ul-20.gif', 'Show/Hide Legend');
hx.addCustomButton('off', @in_ShowLegend, 'checklist_rtl_18dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'What the Colours and Sizes Mean');

hx.show(parentfig);

% ---------------------------------------------------------------- drawing

    function in_drawpanel(k)
        G = Gs{k};
        p = plot(ax(k), G, 'XData', xy(:, 1), 'YData', xy(:, 2));
        p.NodeColor = nodeColor;
        p.MarkerSize = markerSize;
        p.Marker = nodeMarker;
        p.NodeLabelColor = i_labelcolor(hFig, nodeinfo.Highlight);
        p.NodeFontWeight = i_labelweight(nodeinfo.Highlight);
        p.EdgeColor = i_edgecolor(G);
        p.LineWidth = i_edgewidth(G, nodeinfo.EdgeRef(k), edgeWidthMax);
        p.ButtonDownFcn = @(~, ~) in_clickpanel(k);
        pv(k) = p;
        title(ax(k), nodeinfo.PanelTitles(k), 'Interpreter', 'none');
        i_frameaxes(ax(k));
        i_applylabelmode(p);
    end

    function i_applylabelmode(p)
        switch labelMode
            case 1
                p.NodeLabel = cellstr(nodename);
            case 2
                lbl = repmat({''}, nnode, 1);
                lbl(nodeinfo.Highlight) = cellstr(nodename(nodeinfo.Highlight));
                p.NodeLabel = lbl;
            otherwise
                p.NodeLabel = repmat({''}, nnode, 1);
        end
    end

% ------------------------------------------------------------- selection

    function in_clickpanel(k)
        cp = get(ax(k), 'CurrentPoint');
        idx = dsearchn(xy, cp(1, 1:2));
        if idx == selectedNode
            in_ClearSelection([], []);
            return;
        end
        selectedNode = idx;
        for j = 1:npanel
            i_highlightnode(pv(j), Gs{j}, idx, nodeinfo.EdgeRef(j));
        end
        title(t, sprintf('%s  |  selected: %s', figname, nodename(idx)), ...
            'Interpreter', 'none');
    end

    function i_highlightnode(p, G, idx, edgeref)
        [sIdx, tIdx] = findedge(G);
        incident = (sIdx == idx) | (tIdx == idx);
        cc = i_edgecolor(G);
        cc(~incident, :) = repmat([0.85, 0.85, 0.85], sum(~incident), 1);
        lw = i_edgewidth(G, edgeref, edgeWidthMax);
        if isscalar(lw), lw = repmat(lw, G.numedges, 1); end
        if G.numedges > 0
            lw(~incident) = 0.25;
            lw(incident) = max(lw(incident), edgeWidthMax);
            p.EdgeColor = cc;
            p.LineWidth = lw;
        end
        ms = markerSize;
        ms(idx) = ms(idx) + 6;
        p.MarkerSize = ms;
    end

    function in_ClearSelection(~, ~)
        selectedNode = 0;
        for j = 1:npanel
            pv(j).EdgeColor = i_edgecolor(Gs{j});
            pv(j).LineWidth = i_edgewidth(Gs{j}, nodeinfo.EdgeRef(j), edgeWidthMax);
            pv(j).MarkerSize = markerSize;
        end
        title(t, figname, 'Interpreter', 'none');
    end

% --------------------------------------------------------------- buttons

    function in_ChangeFontSize(~, ~)
        for j = 1:npanel
            if pv(j).NodeFontSize >= 20
                pv(j).NodeFontSize = 7;
            else
                pv(j).NodeFontSize = pv(j).NodeFontSize + 1;
            end
        end
    end

    function in_ChangeWeight(~, ~)
        edgeWidthMax = edgeWidthMax + 1;
        if edgeWidthMax > 6, edgeWidthMax = 1; end
        in_ClearSelection([], []);
    end

    function in_ChangeLabels(~, ~)
        labelMode = labelMode + 1;
        if labelMode > 3, labelMode = 1; end
        for j = 1:npanel
            i_applylabelmode(pv(j));
        end
    end

    function in_ChangeLayout(~, ~)
        list = {'force (default)', 'circle', 'subspace', 'layered'};
        [indx, tf] = gui.myListdlg(hFig, list, 'Network Layout', [], false);
        if ~tf, return; end
        methods = {'force', 'circle', 'subspace', 'layered'};
        xy = i_unionlayout(Gs, nodename, methods{indx});
        for j = 1:npanel
            pv(j).XData = xy(:, 1);
            pv(j).YData = xy(:, 2);
        end
    end

    function in_SaveAdj(~, ~)
        if ismcc || isdeployed
            gui.myErrordlg(hFig, 'Export is not available in the standalone application.');
            return;
        end
        labels = cell(1, npanel+1);
        vars = cell(1, npanel+1);
        values = cell(1, npanel+1);
        for j = 1:npanel
            labels{j} = sprintf('Save adjacency matrix of panel %d (%s) to:', ...
                j, nodeinfo.PanelTitles(j));
            vars{j} = sprintf('A%d', j);
            values{j} = full(adjacency(Gs{j}, 'weighted'));
        end
        labels{end} = 'Save gene list to variable named:';
        vars{end} = 'g';
        values{end} = nodename;
        msgfig = export2wsdlg(labels, vars, values);
        uiwait(msgfig);
    end

    function lgd = in_buildlegend()
        % Proxy handles, because a GraphPlot is one object and legend()
        % needs a handle per entry. Only what is actually on screen gets a
        % row: a legend listing a marker the reader cannot find is worse
        % than no legend.
        grey = i_plaincolor();
        accent = gui.i_accentred(hFig);
        [poscol, negcol] = i_edgecolorends();
        anyneg = false;
        anypos = false;
        for j = 1:npanel
            if Gs{j}.numedges > 0 && ...
                    ismember('Weight', Gs{j}.Edges.Properties.VariableNames)
                anyneg = anyneg || any(Gs{j}.Edges.Weight < 0);
                anypos = anypos || any(Gs{j}.Edges.Weight > 0);
            end
        end

        spec = struct('marker', {}, 'color', {}, 'lw', {}, 'label', {});
        if any(nodeinfo.Highlight)
            spec(end+1) = struct('marker', 'o', 'color', accent, 'lw', 0.5, ...
                'label', char(nodeinfo.HighlightName));
        end
        if any(~nodeinfo.Highlight)
            spec(end+1) = struct('marker', 'o', 'color', grey, 'lw', 0.5, ...
                'label', char(nodeinfo.PlainName));
        end
        if any(nodeinfo.IsTF)
            spec(end+1) = struct('marker', 'd', 'color', grey, 'lw', 0.5, ...
                'label', 'transcription factor');
        end
        if anypos
            spec(end+1) = struct('marker', '-', 'color', poscol, 'lw', 2, ...
                'label', 'positive weight');
        end
        if anyneg
            spec(end+1) = struct('marker', '-', 'color', negcol, 'lw', 2, ...
                'label', 'negative weight');
        end

        target = ax(1);
        hold(target, 'on');
        proxies = gobjects(numel(spec), 1);
        for j = 1:numel(spec)
            proxies(j) = plot(target, NaN, NaN, spec(j).marker, ...
                'Color', spec(j).color, ...
                'MarkerFaceColor', spec(j).color, ...
                'MarkerEdgeColor', spec(j).color, ...
                'LineWidth', spec(j).lw, 'MarkerSize', 7);
        end
        lgd = legend(target, proxies, {spec.label}, 'Interpreter', 'none', ...
            'Orientation', 'horizontal', 'FontSize', 8, 'Box', 'off');
        lgd.Layout.Tile = 'south';
    end

    function in_ToggleLegend(~, ~)
        if ~isempty(lgd) && isvalid(lgd)
            delete(lgd);
            lgd = [];
        else
            lgd = in_buildlegend();
        end
    end

    function in_ShowLegend(~, ~)
        msg = sprintf(['Node size and colour: %s (larger and redder = higher).\n', ...
            'Red bold label: seed gene.\n', ...
            'Diamond marker: transcription factor.\n', ...
            'Blue edge: positive weight.  Red edge: negative weight.\n', ...
            'Edge width: |weight| relative to the panel reference.\n\n', ...
            'Click a node to highlight its edges in every panel; ', ...
            'click it again to clear.\n%s'], ...
            nodeinfo.ScoreName, nodeinfo.Legend);
        gui.myHelpdlg(hFig, msg, 'Network Legend');
    end
end

% ------------------------------------------------------- local functions

function xy = i_unionlayout(Gs, nodename, method)
% Lay out the union of all panels: an edge present in any panel pulls its
% endpoints together, so a panel that lost an edge keeps its nodes in place.
if nargin < 3 || isempty(method), method = 'force'; end
nnode = numel(nodename);
U = sparse(nnode, nnode);
for j = 1:numel(Gs)
    U = U + abs(adjacency(Gs{j}, 'weighted'));
end
Gu = graph(max(U, U.'), cellstr(nodename), 'omitselfloops');
fTmp = figure('Visible', 'off');
cleaner = onCleanup(@() close(fTmp));
pTmp = plot(axes(fTmp), Gu);
if strcmp(method, 'force')
    layout(pTmp, 'force', 'Iterations', 300, 'UseGravity', true);
else
    layout(pTmp, method);
end
xy = [pTmp.XData(:), pTmp.YData(:)];
end

function nodeinfo = i_filldefaults(nodeinfo, npanel, nnode)
defaults = struct( ...
    'PanelTitles', "Panel " + string(1:npanel), ...
    'Highlight', false(nnode, 1), ...
    'HighlightName', "highlighted", ...
    'PlainName', "other", ...
    'IsTF', false(nnode, 1), ...
    'EdgeRef', nan(1, npanel), ...
    'Legend', "");
fn = fieldnames(defaults);
for k = 1:numel(fn)
    if ~isfield(nodeinfo, fn{k}) || isempty(nodeinfo.(fn{k}))
        nodeinfo.(fn{k}) = defaults.(fn{k});
    end
end
nodeinfo.PanelTitles = string(nodeinfo.PanelTitles(:))';
nodeinfo.Highlight = logical(nodeinfo.Highlight(:));
nodeinfo.IsTF = logical(nodeinfo.IsTF(:));
nodeinfo.EdgeRef = double(nodeinfo.EdgeRef(:))';
end

function grey = i_plaincolor()
% Shared with the legend so the swatch and the nodes cannot drift apart.
grey = [0.72, 0.72, 0.72];
end

function cc = i_nodecolor(highlight, accent)
% Binary, not a ramp: the fill answers "is this one significant", and node
% SIZE carries the magnitude. Encoding the same score in both left two
% channels saying one thing and no channel saying the other.
cc = repmat(i_plaincolor(), numel(highlight), 1);
cc(highlight, :) = repmat(accent, sum(highlight), 1);
end

function i_frameaxes(ax)
% A frame around each panel, with no numbers on it. Force-layout
% coordinates are arbitrary, so tick values would invite a reading of
% distance the layout does not support - but the frame still separates the
% panels. Colours are left alone so the axes keeps following the figure
% theme, light or dark.
set(ax, 'Box', 'on', 'Layer', 'top', 'XTick', [], 'YTick', []);
end

function cc = i_labelcolor(hFig, isseed)
bkcolor = gui.i_getthemebkgcolor(hFig);
cc = repmat(1-bkcolor, numel(isseed), 1);
cc(isseed, :) = repmat(gui.i_accentred(hFig), sum(isseed), 1);
end

function fw = i_labelweight(isseed)
fw = repmat({'normal'}, numel(isseed), 1);
fw(isseed) = {'bold'};
end

function [pos, neg] = i_edgecolorends()
% Shared with the legend, for the same reason as i_plaincolor.
pos = [0, 0.4470, 0.7410];
neg = [0.8500, 0.3250, 0.0980];
end

function cc = i_edgecolor(G)
[pos, negc] = i_edgecolorends();
n = G.numedges;
cc = repmat(pos, n, 1);
if n > 0 && ismember('Weight', G.Edges.Properties.VariableNames)
    neg = G.Edges.Weight < 0;
    cc(neg, :) = repmat(negc, sum(neg), 1);
end
end

function lw = i_edgewidth(G, ref, lwmax)
n = G.numedges;
if n == 0, lw = 0.5; return; end
if ~ismember('Weight', G.Edges.Properties.VariableNames)
    lw = repmat(0.5*lwmax, n, 1);
    return;
end
a = abs(G.Edges.Weight);
if ~isfinite(ref) || ref <= 0, ref = max(a); end
if ~isfinite(ref) || ref <= 0
    lw = repmat(0.5*lwmax, n, 1);
else
    % A low floor matters more than a high ceiling: these graphs are dense
    % enough that a heavy stroke on every weak edge fills the cluster in
    % and buries the labels.
    lw = 0.25 + (lwmax - 0.25) * min(a/ref, 1);
end
end
