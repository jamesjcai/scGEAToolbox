function xy = i_causalcccdeclutter(xy, side)
%GUI.I_CAUSALCCCDECLUTTER  Enforce a minimum y-gap between same-side nodes,
%so labels don't overlap in GUI.I_NETWORKVISCCC's print view.
%
%   xy = GUI.I_CAUSALCCCDECLUTTER(xy, side)
%
%   xy    - Nx2, as returned by GUI.I_CAUSALCCCGRAPH's plot
%           ([p.XData', p.YData'])
%   side  - Nx1 string, G.Nodes.Side ("senders"/"receivers"); an empty
%           entry is left out
%
%   GUI.I_CAUSALCCCGRAPH's force layout only constrains the ligand/receptor
%   anchors to line up across the two sides (its I_ANCHORALIGN); the
%   intracellular genes around them are free to land wherever the force
%   layout and gravity happened to put them, which on a busy network can
%   leave several nodes within a very short span of y - fine for a
%   GraphPlot's small markers, but GUI.I_NETWORKVISCCC's label-as-marker
%   print view has no marker to fall back on, so overlapping text is
%   unreadable there even though the same xy looks fine as a GraphPlot.
%
%   Standard two-pass minimum-gap enforcement, once per side (senders and
%   receivers do not compete for the same vertical room, since they are
%   drawn in separate x bands): sort by y, push forward, then push
%   backward, so every node keeps at least a minimum gap from its
%   same-side neighbours. Order is preserved exactly; only spacing
%   changes. Rescaled back into [-1,1] afterwards (skipped for a side
%   whose nodes end up fully coincident, which would rescale to NaN),
%   since enforcing gaps can push the total span wider than the original
%   layout's. MINGAP is relative to an even spread of that side's node
%   count across the full range, not a fixed number, so a sparse side is
%   left close to its original (already fine) spacing and only a
%   genuinely crowded side gets pulled apart.
%
%   This does not touch x. A related but separate problem - a label near
%   either side's inner x edge spanning the gap and colliding with one on
%   the OTHER side - is fixed once, upstream, in
%   GUI.I_CAUSALCCCGRAPH>I_BIPARTITELAYOUT's own band widths, since that
%   layout is what every consumer of the graph (this GraphPlot and
%   GUI.I_NETWORKVISCCC alike) draws from; duplicating a second x-fix here
%   would only add a step that no longer does anything.
%
%   USAGE
%     p = gui.i_causalcccgraph(G, figname, []);
%     xy = gui.i_causalcccdeclutter([p.XData', p.YData'], string(G.Nodes.Side));
%     gui.i_networkvisccc(G, xy, true, p.NodeFontSize, []);
%
%   See also GUI.I_CAUSALCCCGRAPH, GUI.I_NETWORKVISCCC.

MINGAP_SCALE = 1.7;

sides = unique(side(side ~= ""), 'stable');
for k = 1:numel(sides)
    idx = find(side == sides(k));
    n = numel(idx);
    if n < 2, continue; end
    mingap = MINGAP_SCALE / (n - 1);

    [ys, ord] = sort(xy(idx, 2));
    for i = 2:n
        if ys(i) - ys(i-1) < mingap
            ys(i) = ys(i-1) + mingap;
        end
    end
    for i = n-1:-1:1
        if ys(i+1) - ys(i) < mingap
            ys(i) = ys(i+1) - mingap;
        end
    end
    if max(ys) - min(ys) > eps
        ys = rescale(ys, -1, 1);
    end

    yout = zeros(n, 1);
    yout(ord) = ys;
    xy(idx, 2) = yout;
end
end
