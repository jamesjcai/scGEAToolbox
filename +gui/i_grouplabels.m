function [groups, labels] = i_grouplabels(values)
%I_GROUPLABELS Groups of a cell grouping, largest first, with their sizes.
%   [GROUPS, LABELS] = gui.i_grouplabels(VALUES) returns the distinct
%   values of the string column VALUES ordered by size, largest first, and
%   LABELS of the form "Beta cells (1229 cells)" for a pick list.
%
%   See also gui.i_pickgrncells, gui.i_cellgroupings.

[gid, groups] = findgroups(values);
counts = accumarray(gid(:), 1);
[counts, order] = sort(counts, 'descend');
groups = groups(order);
groups = groups(:);
labels = groups + " (" + string(counts(:)) + " cells)";
end
