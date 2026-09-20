function h = i_drawgroupmap(ax, h, spec)
%I_DRAWGROUPMAP Draw a gene-by-cell map with its group blocks marked.
%
%   h = gui.i_drawgroupmap(ax, h, spec)
%
%   Inputs:
%     ax   - axes to draw into. Cleared of the previous image and of the
%            previous group separators, nothing else.
%     h    - the image drawn last time, deleted first. GOBJECTS(0) on the
%            first call.
%     spec - a struct describing the map:
%              Y      genes-by-cells matrix, cells already grouped
%              cL     the group labels, in drawing order
%              glist  the gene names, one per row of Y
%              a, b   group tick positions: the ticks go at a-b
%              szgn   how many cells are in each group, in drawing order
%              fliped true to draw the map transposed, genes across
%
%   Returns the new image handle.
%
%   One function because the separators have to be redrawn with the map
%   every time, and nothing in the axes remembers them. IMAGESC goes
%   through NEWPLOT, and an axes left at the default NextPlot='replace'
%   is emptied by it -- image, lines and all. The separators were drawn
%   once at startup, so 'Flip XY' and 'Change normalization method...',
%   which both called IMAGESC and drew no lines of their own, wiped them
%   and never put them back: from the first click on either, the map
%   showed no group boundaries at all for the rest of its life.
%
%   The DELETEs below are therefore belt and braces rather than the
%   mechanism. They are kept so the clearing is written down instead of
%   depending on a NextPlot the caller could reasonably change.
%
%   Both branches state their tick label rotation rather than only the
%   flipped one, which wants its gene names on their side. Setting
%   XTickLabel happens to put the rotation back to 0 by itself, so the
%   upright branch's 0 changes nothing today; it is here so the drawn
%   state is written down rather than inherited.
%
%   See also GUI.I_HEATMAP, GUI.I_GROUPORDERPERM.

arguments
    ax (1,1) matlab.graphics.axis.Axes
    h
    spec (1,1) struct
end

delete(h);
delete(findall(ax, 'Type', 'ConstantLine'));

szc = cumsum(spec.szgn);
ticks = spec.a - spec.b;

if spec.fliped
    h = imagesc(ax, spec.Y.');
    set(ax, 'YTick', ticks);
    set(ax, 'YTickLabel', gui.i_escapeunderscore(spec.cL));
    set(ax, 'XTick', 1:length(spec.glist));
    set(ax, 'XTickLabel', spec.glist);
    set(ax, 'XTickLabelRotation', 90);
    for k = 1:length(szc)
        yline(ax, szc(k)+0.5, 'y-');
    end
else
    h = imagesc(ax, spec.Y);
    set(ax, 'XTick', ticks);
    set(ax, 'XTickLabel', gui.i_escapeunderscore(spec.cL));
    set(ax, 'XTickLabelRotation', 0);
    set(ax, 'YTick', 1:length(spec.glist));
    set(ax, 'YTickLabel', spec.glist);
    for k = 1:length(szc)
        xline(ax, szc(k)+0.5, 'y-');
    end
end

set(ax, 'TickLength', [0, 0]);
box(ax, 'on');

end
