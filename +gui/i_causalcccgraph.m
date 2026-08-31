function p = i_causalcccgraph(G, figname, parentfig)
%GUI.I_CAUSALCCCGRAPH  Render a network built by SC_CAUSALCCCNET.
%
%   p = GUI.I_CAUSALCCCGRAPH(G, figname, parentfig)
%
%   figname is the title. Pass a 2-element string array [title; subtitle]
%   to keep the identifying part (e.g. "microglia -> tumor cell") as a
%   short title with the rest (cell/node/edge counts) as a smaller
%   subtitle underneath, rather than one long title line; a scalar string
%   still works and is shown as a plain title with no subtitle.
%
%   Nodes are coloured by G.Nodes.Color/Group (ligand, sender_gene,
%   receptor, receiver_gene, metadata - the same palette PKG.SC_SCE2CAUSALCCC's
%   state-order file uses). Sender-side nodes (G.Nodes.Side=="senders") are
%   laid out on the left, receiver-side on the right, each by its own force
%   layout so the ligand-receptor bridge edges are not what drags the two
%   sides together - a plain whole-graph force layout tends to pull
%   strongly-linked bridge pairs into one hairball and lose the senders/
%   receivers separation that is the point of this view. Bridge edges
%   (G.Edges.Type=="bridge") are drawn thicker and gold; inferred
%   intracellular edges are thin, coloured by sign where the chosen Method
%   is signed and grey otherwise.
%
%   Nodes can be dragged with the mouse (as in gui.i_singlegraph). The
%   toolbar also has Change Font Size and Change Width of Edges buttons.
%   Layout/directed/cutoff-style controls are deliberately not offered
%   here: those rebuild the graph from a plain adjacency matrix in
%   gui.i_singlegraph, which would silently discard the Side/Type/Sign
%   columns this view's coloring depends on.

if nargin < 3, parentfig = []; end
if nargin < 2, figname = ''; end

hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;
h1 = axes(hFig);

w = 3;
oldidx = 0;

xy = i_bipartitelayout(G, h1);

p = plot(h1, G, 'XData', xy(:, 1), 'YData', xy(:, 2), 'ArrowSize', 6, ...
    'ButtonDownFcn', @startDragFcn);
p.MarkerSize = 6;
% GraphPlot's default Interpreter is 'tex', which reads an underscore as a
% subscript marker - node names routinely have one (i_disambiguate suffixes
% a symbol that appears on both sides with "_senders"/"_receivers", and
% metadata nodes are named "<attr>_senders"/"<attr>_receivers"), so without
% this a name like RELB_senders would render as RELB with a subscript s.
p.Interpreter = 'none';
p.NodeFontSize = 12;   % default GraphPlot NodeFontSize (10) reads small
                       % against this many overlapping node labels

rgb = zeros(G.numnodes, 3);
for i = 1:G.numnodes
    rgb(i, :) = i_hex2rgb(G.Nodes.Color(i));
end
p.NodeColor = rgb;

isBridge = G.Edges.Type == "bridge";
ec = repmat([0.55, 0.55, 0.55], G.numedges, 1);
ec(G.Edges.Sign < 0, :) = repmat([0.85, 0.33, 0.10], sum(G.Edges.Sign < 0), 1);
ec(isBridge, :) = repmat([0.85, 0.65, 0.05], sum(isBridge), 1);
p.EdgeColor = ec;

lw = ones(G.numedges, 1);
inferred = ~isBridge;
if any(inferred) && max(G.Edges.Weight(inferred)) > 0
    lw(inferred) = rescale(G.Edges.Weight(inferred), 1, 4);
end
lw(isBridge) = 3;
p.LineWidth = lw;

i_printlegend(G, isBridge);

hx.addCustomButton('off', @ChangeFontSize, 'noun_font_size_591141.gif', 'Change Font Size of Nodes');
hx.addCustomButton('on',  @ChangeWeight, 'weight_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Width of Edges');

figname = string(figname);
title(h1, figname(1), 'Interpreter', 'none');
if numel(figname) > 1
    subtitle(h1, strjoin(figname(2:end), newline), 'Interpreter', 'none');
end
hx.show(parentfig);


    function ChangeFontSize(~, ~)
        if p.NodeFontSize >= 20
            p.NodeFontSize = 7;
        else
            p.NodeFontSize = p.NodeFontSize + 1;
        end
    end

    function ChangeWeight(~, ~)
        w = w + 1;
        if w > 10, w = 2; end
        p.LineWidth = rescale(p.LineWidth, 1, w);
    end

    % Callback to initiate dragging (mirrors gui.i_singlegraph)
    function startDragFcn(hObj, ~)
        set(hFig, 'WindowButtonMotionFcn', {@draggingFcn, hObj});
        set(hFig, 'WindowButtonUpFcn', @stopDragFcn);
    end

    % Function to drag the point
    function draggingFcn(~, ~, hObj)
        cp = get(gca, 'CurrentPoint');
        yData = get(hObj, 'YData');
        xData = get(hObj, 'XData');
        if oldidx == 0
            idx = dsearchn([xData' yData'], [cp(1,1) cp(1,2)]);
            oldidx = idx;
        else
            idx = oldidx;
        end
        xData(idx) = cp(1,1);
        yData(idx) = cp(1,2);
        set(hObj, 'XData', xData);
        set(hObj, 'YData', yData);
    end

    % Function to stop dragging
    function stopDragFcn(~, ~)
        set(hFig, 'WindowButtonMotionFcn', '');
        set(hFig, 'WindowButtonUpFcn', '');
        oldidx = 0;
    end
end


function xy = i_bipartitelayout(G, h1)
% Force-layout each side on its own, then rescale sender-side x into
% [-1,-0.15] and receiver-side x into [0.15,1], so the bridge edges
% visibly cross a gap instead of being buried inside one hairball.
xy = zeros(G.numnodes, 2);
isSnd = G.Nodes.Side == "senders";
isRcv = G.Nodes.Side == "receivers";

if any(isSnd)
    xy(isSnd, :) = i_sidelayout(subgraph(G, find(isSnd)), h1, -1, -0.15);
end
if any(isRcv)
    xy(isRcv, :) = i_sidelayout(subgraph(G, find(isRcv)), h1, 0.15, 1);
end
end


function xy = i_sidelayout(Gside, h1, xlo, xhi)
p0 = plot(h1, Gside);
% UseGravity matters here specifically because a side's intracellular
% network is often not fully connected (TopKPerNode sparsification and
% weak MI/correlation edges both leave isolated nodes or split it into
% several small components). Without it, layout('force') snaps any graph
% with more than one component onto a rigid grid, one cell per component -
% which reads as a different "style" of layout than a fully connected
% side gets, even though both call layout with the same method. Gravity
% pulls every node toward a shared centroid regardless of connectivity, so
% both sides always get an organic force spread instead of one of them
% occasionally turning into a grid.
layout(p0, 'force', 'UseGravity', true);
x = rescale(p0.XData(:), xlo, xhi);
y = rescale(p0.YData(:), -1, 1);
delete(p0);
xy = [x, y];
end


function rgb = i_hex2rgb(hexstr)
hexstr = char(erase(hexstr, "#"));
rgb = [hex2dec(hexstr(1:2)), hex2dec(hexstr(3:4)), hex2dec(hexstr(5:6))] / 255;
end


function i_printlegend(G, isBridge)
% Text legend on the command line for the node/edge colors this function
% just drew - there is no on-figure legend, since a Group/Type/Sign
% legend box would just be one more thing crowding an already dense plot.
% Only the rows that are actually present in G are printed, so a network
% with no negative edges (an unsigned Method, or one that simply has none)
% is not told about a color it does not use.
nodeInfo = struct( ...
    'ligand',        {{"dark green (#2e6e5e)",  "selected ligand, sender side"}}, ...
    'sender_gene',   {{"light green (#a9dcc9)", "other sender-side intracellular gene"}}, ...
    'receptor',      {{"brown (#a05a2c)",       "selected receptor, receiver side"}}, ...
    'receiver_gene', {{"cream (#f3e6d8)",       "other receiver-side intracellular gene"}}, ...
    'metadata',      {{"grey (#8b938c)",        "categorical metadata pivot"}});

fprintf('\n[causalCCC network legend]\n');
fprintf('Node color (by G.Nodes.Group):\n');
present = unique(G.Nodes.Group, 'stable');
for grp = present(:)'
    key = char(grp);
    if isfield(nodeInfo, key)
        fprintf('  %-16s %-24s - %s\n', grp, nodeInfo.(key){1}, nodeInfo.(key){2});
    else
        fprintf('  %-16s %-24s\n', grp, '(unrecognized group)');
    end
end

fprintf('Edge color:\n');
if any(isBridge)
    fprintf('  %-16s %-24s - %s\n', 'bridge', 'gold, thick', ...
        'known ligand-receptor link (G.Edges.Type=="bridge"), not inferred');
end
inferred = ~isBridge;
if any(inferred & G.Edges.Sign >= 0)
    fprintf('  %-16s %-24s - %s\n', 'inferred', 'grey, thin', ...
        'intracellular dependency edge, non-negative or from an unsigned Method');
end
if any(inferred & G.Edges.Sign < 0)
    fprintf('  %-16s %-24s - %s\n', 'inferred', 'orange, thin', ...
        'intracellular dependency edge, negative (signed Method only)');
end
fprintf('\n');
end
