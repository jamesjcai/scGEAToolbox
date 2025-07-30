function p = i_violinplot_groupordered(d, c, sorted_c, colorit)
% UNUSED function will be removed.

import pkg.Violin
import pkg.violinplot

if nargin < 4, colorit = false; end
if nargin < 3, sorted_c = []; end
if ~isstring(c), c = string(c); end

c = strrep(c, '_', ' ');
[~, cL] = findgroups(string(c));
if isempty(sorted_c)
    % [~, sortidx] = sort(grpstats(d, c, @median), 'descend');
    [~, sortidx] = sort(splitapply(@median, d, c), 'descend');
    sorted_c = cL(sortidx);
else
    sorted_c = strrep(sorted_c, '_', ' ');
end

if ~iscell(sorted_c)
    sorted_c = cellstr(sorted_c);
end
if colorit
    p = violinplot(d, c, 'GroupOrder', sorted_c, ...
        'ShowData', false, ...
        'EdgeColor', [0, 0, 0]);
else
    p = violinplot(d, c, 'GroupOrder', sorted_c, ...
        'ShowData', false, 'ViolinColor', [1, 1, 1], ...
        'EdgeColor', [0, 0, 0]);
end
xtickangle(-45);
box on
ylabel('Differentiation Potency');
grid on
