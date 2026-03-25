function i_bindi_violinplot_base(d, c, colorit, grouporder)
import gui.Violin
import gui.i_violinplot_base

if ~isstring(c)
    c = string(c);
end
if nargin < 4, grouporder = []; end
if nargin < 3 || isempty(colorit), colorit = false; end
if issparse(d), d = full(d); end

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
% xtickangle(-45);
box on
grid on
end
