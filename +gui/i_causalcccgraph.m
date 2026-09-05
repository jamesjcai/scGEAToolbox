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
%   state-order file uses). Ligand and receptor nodes are drawn larger than
%   the intracellular genes, since they are the anchors the whole view is
%   built around. Sender-side nodes (G.Nodes.Side=="senders") are laid out
%   on the left, receiver-side on the right, each by its own force layout so
%   the ligand-receptor bridge edges are not what drags the two sides
%   together - a plain whole-graph force layout tends to pull strongly-linked
%   bridge pairs into one hairball and lose the senders/receivers separation
%   that is the point of this view. Each side is then rotated to line its
%   bridge endpoints up with the other side (see I_ALIGNSIDES), and finally
%   passed through GUI.I_CAUSALCCCDECLUTTER, since even this GraphPlot's own
%   node names are text that can land close enough in y to overlap - the
%   force layout and I_ANCHORALIGN's partial nudge do not guarantee any
%   minimum spacing on their own. Bridge edges (G.Edges.Type=="bridge") are
%   drawn thicker and gold; inferred intracellular edges are thin, coloured
%   by sign where the chosen Method is signed and grey otherwise.
%
%   The axes is drawn as a plain framed panel with no tick marks: layout
%   coordinates carry no meaning, but the axes has to stay visible for its
%   background to be painted, so the panel follows the figure theme (white
%   under the light one) instead of showing the figure colour through.
%
%   The two Network Vis buttons hand the current layout - including
%   anything dragged since this figure was drawn, since the handoff reads
%   the coordinates off the plot rather than recomputing them - through
%   GUI.I_CAUSALCCCDECLUTTER again (a no-op unless something was dragged
%   into a new collision) to GUI.I_NETWORKVISCCC, whose label-centred
%   drawing prints better than a GraphPlot. That view styles itself from
%   the same I_CCCSTYLE, so the print reads the same way as this screen.
%
%   Nodes can be dragged with the mouse (as in gui.i_singlegraph). Hovering
%   a node gives its group and side, which is how to read a label that has
%   been overplotted by its neighbours. The toolbar has Change Font Size,
%   Change Width of Edges, and Show/Hide Legend buttons.
%   Layout/directed/cutoff-style controls are deliberately not offered
%   here: those rebuild the graph from a plain adjacency matrix in
%   gui.i_singlegraph, which would silently discard the Side/Type/Sign
%   columns this view's coloring depends on.
%
%   See also SC_CAUSALCCCNET, GUI.E_SCE2CAUSALCCCVIEW.

if nargin < 3, parentfig = []; end
if nargin < 2, figname = ''; end

% forceclassic is passed explicitly, not because the default differs - it is
% already true - but because it is load-bearing here and the class help says
% the opposite ("uifigure parent -> uifigure child"). The node dragging below
% needs WindowButtonMotionFcn and the axes CurrentPoint, which do not work on
% a uifigure, so if that default ever changes this view must not follow it.
hx = gui.myFigure(parentfig, true);
hFig = hx.FigHandle;
% myFigure already made the axes; calling axes(hFig) again would leave a
% second, empty one on top, whose ruler is then drawn over this one's -
% which is where the doubled, overlapping tick labels came from.
h1 = hx.AxHandle;

oldidx = 0;
lgd = [];

xy = i_bipartitelayout(G, h1);
% This GraphPlot's own node names are text too, not just markers - a
% same-side pair that lands within a fraction of a unit of each other in y
% (I_ANCHORALIGN's NUDGE is partial by design, so even an anchor is not
% pinned exactly) prints as overlapping, unreadable text right here, not
% only in GUI.I_NETWORKVISCCC's handoff. Same fix, same place it is needed.
xy = gui.i_causalcccdeclutter(xy, string(G.Nodes.Side));

p = plot(h1, G, 'XData', xy(:, 1), 'YData', xy(:, 2), 'ArrowSize', 6, ...
    'ButtonDownFcn', @startDragFcn);
% GraphPlot's default Interpreter is 'tex', which reads an underscore as a
% subscript marker - node names routinely have one (i_disambiguate suffixes
% a symbol that appears on both sides with "_senders"/"_receivers", and
% metadata nodes are named "<attr>_senders"/"<attr>_receivers"), so without
% this a name like RELB_senders would render as RELB with a subscript s.
p.Interpreter = 'none';
p.NodeFontSize = 12;   % default GraphPlot NodeFontSize (10) reads small
                       % against this many overlapping node labels

% Every colour and width comes from i_cccstyle, which reads them off the
% graph's own columns. gui.i_networkvisccc uses the same helper, so the print
% view and this one cannot drift apart on what a colour means.
%
% The ligand and receptor anchors carry the interaction and are drawn larger
% than the intracellular genes, which are context.
st = i_cccstyle(G);
isBridge = st.IsBridge;
p.MarkerSize = st.NodeSize;
p.NodeColor = st.NodeColor;
p.EdgeColor = st.EdgeColor;

% baseLW is kept so ChangeWeight can rescale from the original widths every
% time. Rescaling p.LineWidth in place compounds - each press stretches an
% already-stretched range, and after a few presses the widths no longer
% reflect the weights at all.
baseLW = st.EdgeWidth;
p.LineWidth = baseLW;
w = 4;

i_adddatatips(p, G);
i_sideheaders(h1, G);

% The ruler is what is worth removing here, not the axes. A layout
% coordinate means nothing in a network, so the ticks and their numbers go -
% but the axes itself stays visible, because an invisible axes does not
% paint its own background and the figure's grey then shows through the plot
% area. Leaving the axes to paint keeps the panel white under the light
% theme, and correctly dark under a dark one, without hard-coding either.
box(h1, 'on');
h1.XTick = [];
h1.YTick = [];
i_padlimits(h1, xy);

i_printlegend(G, isBridge);

hx.addCustomButton('off', @ChangeFontSize, 'noun_font_size_591141.gif', 'Change Font Size of Nodes');
hx.addCustomButton('on',  @ChangeWeight, 'weight_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Width of Edges');
hx.addCustomButton('off', @ToggleLegend, 'icon-fa-list-ul-20.gif', 'Show/Hide Legend');
hx.addCustomButton('on',  @in_networkvis_linear, "linear.jpg", "Stright Network Vis");
hx.addCustomButton('off', @in_networkvis_curvy, "curve-array.jpg", "Curvy Network Vis");

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
        w = w + 2;
        if w > 12, w = 2; end
        p.LineWidth = rescale(baseLW, 1, w);
    end

    function ToggleLegend(~, ~)
        % Off by default: a Group/Type legend box is one more thing crowding
        % an already dense plot, and the same information goes to the command
        % line. On demand it beats going back to read that text.
        if ~isempty(lgd) && isvalid(lgd)
            delete(lgd);
            lgd = [];
            return
        end
        lgd = i_onfigurelegend(h1, G, isBridge);
    end

    function in_networkvis_curvy(~, ~)
        i_handoff(true);
    end

    function in_networkvis_linear(~, ~)
        i_handoff(false);
    end

    function i_handoff(curved)
        % The print view, handed the layout this figure is currently
        % showing - reading p.XData/p.YData rather than recomputing. xy was
        % already run through gui.i_causalcccdeclutter once, above, before
        % this GraphPlot was even drawn; re-running it here is only for
        % whatever the user has dragged since - decluttering and "keep my
        % dragged position exactly" are contradictory goals, and this view
        % is print-oriented output, not a further editing surface, so
        % legibility wins there. Drag on THIS figure's GraphPlot still
        % works and is unaffected; only the handoff copy is moved. It
        % styles itself from G through the same i_cccstyle, so there is
        % nothing else to pass.
        fw = gui.myWaitbar(hFig);
        xyPrint = gui.i_causalcccdeclutter([p.XData' p.YData'], string(G.Nodes.Side));
        gui.i_networkvisccc(G, xyPrint, curved, p.NodeFontSize, hFig);
        gui.myWaitbar(hFig, fw);
    end

    % Callback to initiate dragging (mirrors gui.i_singlegraph)
    function startDragFcn(hObj, ~)
        set(hFig, 'WindowButtonMotionFcn', {@draggingFcn, hObj});
        set(hFig, 'WindowButtonUpFcn', @stopDragFcn);
    end

    % Function to drag the point
    function draggingFcn(~, ~, hObj)
        % h1, not gca: gca follows whatever figure the user last clicked,
        % so a second open figure would send these coordinates to the wrong
        % axes - and on a uifigure gca does not resolve to this axes at all.
        cp = get(h1, 'CurrentPoint');
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
% Force-layout each side on its own, orient the two to face each other,
% then rescale sender-side x into [SND_LO,SND_HI] and receiver-side x into
% [RCV_LO,RCV_HI], so the bridge edges visibly cross a gap instead of being
% buried inside one hairball.
%
% Both bands are 1.5 wide with a 0.7 gap between their inner edges, not the
% original 0.85-wide/0.3-gap. I_UNIT (in I_SIDELAYOUT) already normalizes
% each side's own x AND y independently into [-1,1], and a side's raw
% force layout is close to circular before that - measured on a real
% network, x range 3.958 vs y range 3.968, essentially 1:1 - so squeezing
% that circle's x down to a 0.85-wide band while its y keeps the full
% [-1,1] (height 2) drew every side at a forced 2.35-tall:1-wide ratio
% regardless of what the layout actually looked like, which is what made
% every network read as two narrow vertical columns rather than two
% roughly circular clusters. 1.5 wide brings that down to a much less
% distorted 1.33:1. The 0.7 gap (up from 0.3) is independently needed
% because GUI.I_NETWORKVISCCC's label-as-marker print view has no marker
% to fall back on the way this GraphPlot's small circles do: a long,
% centre-aligned label near either band's inner edge can span a narrow
% gap outright and collide with a similar-height label on the other side
% - seen directly in one network, where a sender at the old x=-0.197 and
% a receiver at the old x=0.150 printed on top of each other. Both fixes
% land here, in the one function that lays out every consumer of this
% graph (this GraphPlot, and GUI.I_NETWORKVISCCC via GUI.I_HANDOFF and
% GUI.I_CAUSALCCCDECLUTTER), rather than only in the print-view path.
SND_LO = -1.85; SND_HI = -0.35;
RCV_LO = 0.35;  RCV_HI = 1.85;

xy = zeros(G.numnodes, 2);
sndIdx = find(G.Nodes.Side == "senders");
rcvIdx = find(G.Nodes.Side == "receivers");

A = i_sidelayout(G, sndIdx, h1);
B = i_sidelayout(G, rcvIdx, h1);
[ia, ib] = i_bridgepairs(G, sndIdx, rcvIdx);
[A, B] = i_alignsides(A, B, ia, ib);
[A, B] = i_anchoralign(G, sndIdx, rcvIdx, A, B, ia, ib);

if ~isempty(sndIdx)
    xy(sndIdx, :) = [i_band(A(:, 1), SND_LO, SND_HI), A(:, 2)];
end
if ~isempty(rcvIdx)
    xy(rcvIdx, :) = [i_band(B(:, 1), RCV_LO, RCV_HI), B(:, 2)];
end
end


function xy = i_sidelayout(G, idx, h1)
if isempty(idx)
    xy = zeros(0, 2);
    return
end
Gside = subgraph(G, idx);
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
xy = i_unit([p0.XData(:), p0.YData(:)]);
delete(p0);
end


function [ia, ib] = i_bridgepairs(G, sndIdx, rcvIdx)
% Bridge endpoints as positions within each side's own node list. A bridge
% may be stored either way round, so both orientations are checked.
ia = zeros(0, 1);
ib = zeros(0, 1);
e = G.Edges(G.Edges.Type == "bridge", :);
if isempty(e) || isempty(sndIdx) || isempty(rcvIdx), return; end

sndPos = zeros(G.numnodes, 1);
sndPos(sndIdx) = 1:numel(sndIdx);
rcvPos = zeros(G.numnodes, 1);
rcvPos(rcvIdx) = 1:numel(rcvIdx);

s = findnode(G, e.EndNodes(:, 1));
t = findnode(G, e.EndNodes(:, 2));
fwd = sndPos(s) > 0 & rcvPos(t) > 0;
rev = sndPos(t) > 0 & rcvPos(s) > 0;
ia = [sndPos(s(fwd)); sndPos(t(rev))];
ib = [rcvPos(t(fwd)); rcvPos(s(rev))];
end


function [A, B] = i_alignsides(A, B, ia, ib)
% A force layout has no preferred orientation - rotating or mirroring one is
% just as valid a drawing of the same graph. So each side is free to turn,
% and the angle worth choosing is the one that puts the two ends of each
% ligand-receptor bridge at a similar height. Without this the bridges leave
% at whatever y the force layout happened to land on and cross each other in
% the middle, which is exactly the part of the figure a reader is trying to
% follow.
%
% The two sides constrain each other, so this alternates: hold one still and
% turn the other, then swap. Three rounds is enough to settle in practice.
if isempty(ia) || isempty(A) || isempty(B), return; end

angles = linspace(0, 2*pi, 49);
angles(end) = [];
for iter = 1:3
    B = i_bestrotation(B, A(ia, 2), ib, angles);
    A = i_bestrotation(A, B(ib, 2), ia, angles);
end
end


function [A, B] = i_anchoralign(G, sndIdx, rcvIdx, A, B, ia, ib)
% Rotating a side helps, but it can only move every node together, and the
% ligands sit wherever their own intracellular network put them. So after
% the rotation the anchors themselves are nudged in y into the order that
% makes their bridges run flat: each anchor is ranked by the mean height of
% the partners it bridges to, and pulled toward an evenly spaced column in
% that order.
%
% Measured on the Macrophages -> Beta cells example over 8 layout seeds,
% this takes the mean bridge |dy| from 0.458 (rotation alone) to 0.102, at
% the cost of about 10% longer intracellular edges - the anchors are pulled
% off their force positions, which is the thing being traded away. NUDGE is
% partial for that reason: it buys the whole alignment gain while leaving
% the anchors nearer where their own neighbours placed them.
%
% Only y moves. x stays as the force layout and the side band set it, so an
% anchor is still drawn inside its own side's cluster.
NUDGE = 0.75;
ROUNDS = 6;   % the two sides define each other's targets, so alternate

if isempty(ia) || isempty(A) || isempty(B), return; end
isLig = G.Nodes.Group(sndIdx) == "ligand";
isRec = G.Nodes.Group(rcvIdx) == "receptor";
for r = 1:ROUNDS
    A = i_nudgeanchors(A, isLig, ia, ib, B, NUDGE);
    B = i_nudgeanchors(B, isRec, ib, ia, A, NUDGE);
end
end


function P = i_nudgeanchors(P, isAnchor, myIdx, othIdx, Q, alpha)
anchors = find(isAnchor);
bary = nan(numel(anchors), 1);
for k = 1:numel(anchors)
    m = myIdx == anchors(k);
    if any(m), bary(k) = mean(Q(othIdx(m), 2)); end
end
% An anchor with no bridge on this side has no opinion about where it
% should go, so it keeps the height the force layout gave it.
ok = ~isnan(bary);
if sum(ok) < 2, return; end
present = anchors(ok);
[~, ord] = sort(bary(ok));
target = linspace(-1, 1, numel(present))';
P(present(ord), 2) = (1 - alpha) * P(present(ord), 2) + alpha * target;
end


function P = i_bestrotation(P, targetY, idx, angles)
% The rotation and reflection of P whose bridge endpoints sit closest in y
% to targetY, judged on the coordinates as they will finally be drawn.
if size(P, 1) < 2, return; end
c = mean(P, 1);
best = inf;
bestP = P;
for refl = [1, -1]
    Q0 = (P - c) .* [1, refl];
    for th = angles
        R = [cos(th), -sin(th); sin(th), cos(th)];
        Q = i_unit(Q0 * R');
        cost = sum((Q(idx, 2) - targetY) .^ 2);
        if cost < best
            best = cost;
            bestP = Q;
        end
    end
end
P = bestP;
end


function P = i_unit(P)
% Both axes into [-1,1] independently, matching how the result is finally
% drawn (x is squeezed into a side band afterwards, y kept as is).
P = [i_band(P(:, 1), -1, 1), i_band(P(:, 2), -1, 1)];
end


function v = i_band(v, lo, hi)
% rescale() returns NaN when every value is identical, which happens for a
% one-node side and for a component that the force layout stacked exactly;
% those belong in the middle of the band, not at NaN.
if isempty(v)
    return
end
if max(v) - min(v) < eps
    v = repmat((lo + hi) / 2, size(v));
else
    v = rescale(v, lo, hi);
end
end


function i_padlimits(h1, xy)
% axis equal fits the data exactly, so labels drawn to the right of the
% rightmost node run outside the axes. Node labels extend to the right, so
% the right margin needs to be the generous one.
if isempty(xy), return; end
% max-min rather than range(), which lives in Statistics and Machine
% Learning Toolbox - not worth a toolbox dependency inside a plot helper.
xr = max(max(xy(:, 1)) - min(xy(:, 1)), eps);
yr = max(max(xy(:, 2)) - min(xy(:, 2)), eps);
xlim(h1, [min(xy(:, 1)) - 0.10 * xr, max(xy(:, 1)) + 0.30 * xr]);
ylim(h1, [min(xy(:, 2)) - 0.08 * yr, max(xy(:, 2)) + 0.12 * yr]);
end


function i_sideheaders(h1, G)
% Which half is which is the first thing a reader needs and the one thing
% the drawing itself does not say. The two x positions are the midpoints
% of I_BIPARTITELAYOUT's SND_LO/SND_HI and RCV_LO/RCV_HI bands - update
% both together if either changes.
y = 1.10;
if any(G.Nodes.Side == "senders")
    text(h1, -1.10, y, 'senders', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 11, 'Color', [0.18, 0.43, 0.37]);
end
if any(G.Nodes.Side == "receivers")
    text(h1, 1.10, y, 'receivers', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 11, 'Color', [0.63, 0.35, 0.17]);
end
end


function i_adddatatips(p, G)
% With this many labels, neighbours overplot each other; hovering is how a
% reader recovers a name that has been buried, so the tip carries the group
% and side too. Wrapped because GraphPlot datatip support has changed across
% releases and a missing tip is not worth failing the whole plot over.
try
    p.DataTipTemplate.DataTipRows = [ ...
        dataTipTextRow('', G.Nodes.Name), ...
        dataTipTextRow('group', G.Nodes.Group), ...
        dataTipTextRow('side', G.Nodes.Side)];
catch
    % older release, or GraphPlot without a DataTipTemplate: skip the tips
end
end


function lgd = i_onfigurelegend(h1, G, ~)
% Proxy line objects, because a GraphPlot is one object and legend() cannot
% pick individual nodes or edges out of it. The rows come from i_cccstyle,
% the same place the colours do.
st = i_cccstyle(G);
info = st.Legend;
hold(h1, 'on');
h = gobjects(numel(info), 1);
for k = 1:numel(info)
    h(k) = plot(h1, NaN, NaN, info(k).marker, ...
        'MarkerFaceColor', info(k).color, 'MarkerEdgeColor', info(k).color, ...
        'Color', info(k).color, 'LineWidth', info(k).lw, 'MarkerSize', 7);
end
hold(h1, 'off');
lgd = legend(h1, h, {info.label}, 'Location', 'southoutside', ...
    'Interpreter', 'none', 'NumColumns', 2, 'FontSize', 8);
lgd.Box = 'off';
end


function i_printlegend(G, isBridge)
% Text legend on the command line for the node/edge colors this function
% just drew. The figure has a Show/Hide Legend button for the same thing;
% this stays because it is also what a scripted, headless call gets.
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
