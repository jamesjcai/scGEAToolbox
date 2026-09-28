function i_bindviolinplot(d, c, colorit, grouporder, ax)
import gui.Violin
import gui.i_violinplot_base

if ~isstring(c)
    c = string(c);
end
if nargin < 5, ax = []; end
if nargin < 4, grouporder = []; end
if nargin < 3 || isempty(colorit), colorit = false; end
if issparse(d), d = full(d); end

% gui.Violin draws into the current axes, so make AX current first. This sets
% the current figure without raising it, unlike axes(AX). A button callback
% cannot rely on gca: after a dialog closes it can be another window's.
if isempty(ax)
    ax = gca;
else
    fig = ancestor(ax, 'figure');
    set(groot, 'CurrentFigure', fig);
    fig.CurrentAxes = ax;
end

if ~colorit
    if isempty(grouporder)
        i_violinplot_base(d, c, ...
            'ShowData', false, 'ViolinColor', [1, 1, 1], ...
            'EdgeColor', [0, 0, 0]);
    else
        if ~iscell(grouporder), grouporder = cellstr(grouporder); end
        i_violinplot_base(d, c, ...
            'ShowData', false, 'ViolinColor', [1, 1, 1], ...
            'EdgeColor', [0, 0, 0], 'GroupOrder', grouporder);
    end
else
    if isempty(grouporder)
        i_violinplot_base(d, c, 'ShowData', false, 'EdgeColor', [0, 0, 0]);
    else
        if ~iscell(grouporder), grouporder = cellstr(grouporder); end
        i_violinplot_base(d, c, 'ShowData', false, 'EdgeColor', [0, 0, 0], ...
            'GroupOrder', grouporder);
    end
end
box(ax, 'on');
grid(ax, 'on');
end
