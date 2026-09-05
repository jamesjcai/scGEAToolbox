function i_networkvisccc(G, xy, curved, fontsize, parentfig)
%GUI.I_NETWORKVISCCC  Print view of a causalCCC network.
%
%   GUI.I_NETWORKVISCCC(G, xy, curved, fontsize, parentfig)
%
%   The same label-centred drawing as GUI.I_NETWORKVIS - node names as the
%   visual, no markers to occlude them, which prints better than a
%   GraphPlot - but styled for a network from SC_CAUSALCCCNET. Colours and
%   widths come from the graph itself via I_CCCSTYLE, so there is nothing to
%   pass and nothing to keep in step by hand:
%
%     - labels carry their group's colour, darkened enough to stay readable
%       (the palette has a cream and a pale green that are fine as filled
%       markers and invisible as letters on white)
%     - known ligand-receptor bridges are gold and thick, always straight,
%       and drawn last so they are not buried under the intracellular
%       hairlines
%     - each label gets an opaque backing, so it reads over the edges
%       beneath rather than having them show through the letters
%     - a legend, since a printed figure has to explain its own colours
%
%   CURVED bends the intracellular edges into quadratic curves; false draws
%   everything straight. Ligand-receptor bridges are straight either way -
%   they are what a reader follows across the gap between the two sides, and
%   a curve there wanders into the opposite cluster. It also leaves the two
%   kinds of edge distinguishable by shape, not colour alone.
%
%   Edge widths are deliberately bridge-or-not only, not scaled by weight:
%   as real lines the weight-scaled widths become a mat dense enough to bury
%   the labels.
%
%   GUI.I_NETWORKVIS stays the general-purpose version and is unaffected by
%   any of this.
%
%   See also GUI.I_NETWORKVIS, GUI.I_CAUSALCCCGRAPH, SC_CAUSALCCCNET.

if nargin < 5, parentfig = []; end
if nargin < 4, fontsize = 15; end
if nargin < 3, curved = true; end

st = i_cccstyle(G);

fx = gui.myFigure(parentfig);
ax = fx.AxHandle;
lgd = [];

% Each drawn pair mapped back to its row in G.Edges: the drawing order comes
% from find() on the adjacency, which is not G.Edges order, so without this
% the styling would land on the wrong edges.
[si, sj] = find(adjacency(G));
[~, ord] = sort(max(si, sj));
si = si(ord);
sj = sj(ord);
eidx = findedge(G, si, sj);

ecol = st.EdgeColor(eidx, :);
ewid = st.PrintWidth(eidx);
ebr = st.IsBridge(eidx);
% Thin first, so the gold bridges end up on top of the grey mesh.
[~, zord] = sort(ewid, 'ascend');
si = si(zord);
sj = sj(zord);
ecol = ecol(zord, :);
ewid = ewid(zord);
ebr = ebr(zord);

i_drawedges();

hold(ax, "on");
scatter(ax, xy(:,1), xy(:,2), st.NodeSize * 4, st.NodeColor, 'filled', ...
    'MarkerEdgeColor', 'none');

i_addtext();
% XColor/YColor 'none' rather than axis off: an invisible axes does not paint
% its own background, so the figure colour would show through the drawing
% area. This hides the ruler and the box while leaving the axes to paint,
% keeping the sheet white under the light theme.
set(ax, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');

lgd = i_legend();
fx.addCustomButton('off', @i_togglelegend, 'icon-fa-list-ul-20.gif', 'Show/Hide Legend');
fx.addCustomButton('off', {@i_rotatetext, true}, "fillet3d.jpg", "Rotate Text");
fx.addCustomButton('off', {@i_rotatetext, false}, "rotation.gif", "Rotate Text");
fx.show(parentfig);
vangle = 0;


    function i_drawedges()
        M = zeros(0, 4);
        hold(ax, "on");
        for e = 1:numel(si)
            p1 = xy(si(e), :);
            p2 = xy(sj(e), :);
            if curved
                % An undirected graph's adjacency lists every edge twice and
                % the second copy would only overprint the first.
                if ismember([p1, p2], M, 'rows'), continue; end
                M = [M; p1, p2; p2, p1]; %#ok<AGROW>
            end
            % A ligand-receptor link is always straight, even in curved mode.
            % It is the one edge a reader follows end to end, and it is the
            % one that crosses the gap between the two sides - a curve there
            % wanders into the opposite cluster and stops being traceable,
            % while the anchor alignment has already made these near
            % horizontal. The intracellular edges still curve, which also
            % makes the two kinds tell themselves apart by shape and not
            % only by colour.
            if curved && ~ebr(e)
                [a, b] = quadraticcurveto(p1, p2);
                plot(ax, a, b, '-', 'color', ecol(e, :), 'linewidth', ewid(e));
            else
                plot(ax, [p1(1), p2(1)], [p1(2), p2(2)], '-', ...
                    'color', ecol(e, :), 'linewidth', ewid(e));
            end
        end
    end

    function i_addtext(rotation)
        if nargin < 1, rotation = 0; end
        % ax, not gcf: gcf is whatever figure the user last clicked, so this
        % could strip the labels off an unrelated figure instead of this one.
        delete(findall(ax, 'Type', 'text'));
        bk = gui.i_getthemebkgcolor(fx.FigHandle);
        for k = 1:G.numnodes
            wx = i_measure(G.Nodes.Name{k});
            text(ax, xy(k,1) - floor(wx/2), xy(k,2), G.Nodes.Name{k}, ...
                'FontSize', fontsize, ...
                'Color', i_readable(st.NodeColor(k, :), bk), ...
                'BackgroundColor', ax.Color, ...
                'FontWeight', 'normal', ...
                'Interpreter', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Margin', 0.2, ...
                'Rotation', rotation);
        end
    end

    function i_rotatetext(~, ~, increase)
        if increase
            vangle = vangle + 5;
        else
            vangle = vangle - 5;
        end
        i_addtext(vangle);
    end

    function wx = i_measure(txt)
        % Interpreter 'none' here too: a TeX-rendered string measures
        % narrower than the literal one, which would throw the centring off.
        h = text(ax, 0, 0, txt, 'FontSize', fontsize, ...
            'FontWeight', 'normal', 'Interpreter', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        ext = get(h, 'Extent');
        delete(h);
        wx = ext(3) / 3;
    end

    function h = i_legend()
        % Proxy lines, because the drawing is hundreds of separate objects
        % and legend() needs one handle per entry.
        info = st.Legend;
        hold(ax, 'on');
        proxies = gobjects(numel(info), 1);
        for k = 1:numel(info)
            proxies(k) = plot(ax, NaN, NaN, info(k).marker, ...
                'MarkerFaceColor', info(k).color, ...
                'MarkerEdgeColor', info(k).color, ...
                'Color', info(k).color, 'LineWidth', info(k).lw, ...
                'MarkerSize', 7);
        end
        h = legend(ax, proxies, {info.label}, 'Location', 'southoutside', ...
            'Interpreter', 'none', 'NumColumns', 2, 'FontSize', 8);
        h.Box = 'off';
    end

    function i_togglelegend(~, ~)
        % On by default here, unlike the interactive view: a figure on paper
        % has to explain its own colours, with no toolbar to fall back on.
        if ~isempty(lgd) && isvalid(lgd)
            delete(lgd);
            lgd = [];
            return
        end
        lgd = i_legend();
    end
end


function c = i_readable(c, bkcolor)
% Keep the hue that identifies a node's group, but move it far enough from
% the background to be legible as text. The causalCCC palette has a cream
% and a pale green that read fine as a filled dot and vanish as letters.
lum = @(x) 0.2126*x(1) + 0.7152*x(2) + 0.0722*x(3);
limit = 0.45;
if lum(bkcolor) > 0.5
    L = lum(c);
    if L > limit
        c = c * (limit / max(L, eps));
    end
else
    L = lum(c);
    if L < 1 - limit
        c = 1 - (1 - c) * (limit / max(1 - L, eps));
    end
end
c = min(max(c, 0), 1);
end
