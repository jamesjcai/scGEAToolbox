function [names, counts] = i_ordergrouplevels(thisc, sortby)
%I_ORDERGROUPLEVELS List the levels of a grouping variable in a given order.
%
%   [names, counts] = gui.i_ordergrouplevels(thisc, sortby)
%
%   thisc  - per-cell group label vector (string, cellstr, or numeric).
%   sortby - "natural" (default), "count", or "none":
%              "natural" alphabetic, with embedded numbers read as numbers,
%                        so Cluster10 sorts after Cluster2 rather than before
%              "count"   most cells first; equal counts keep natural order
%              "none"    the order the levels first appear in the data
%
%   names  - the level names, as a string column, in that order.
%   counts - how many cells carry each, matching NAMES element for element.
%
%   Split out of GUI.I_SELECTGROUPSUBSET so the ordering can be tested
%   without putting a dialog on screen.
%
%   See also GUI.I_SELECTGROUPSUBSET.

if nargin < 2 || isempty(sortby), sortby = "natural"; end
sortby = string(sortby);

thisc = string(thisc);
[ci, cLi] = findgroups(thisc(:));
counts = accumarray(ci(:), 1, [numel(cLi), 1]);
cLi = cLi(:);

switch sortby
    case "count"
        % SORT is stable, so equal counts stay in the alphabetical order
        % FINDGROUPS returned them in.
        [~, ord] = sort(counts, 'descend');
    case "none"
        [~, ia] = unique(thisc(:), 'stable');
        [~, ord] = ismember(thisc(ia), cLi);
    case {"natural", ""}
        [~, ord] = ismember(natsort(cLi), cLi);
    otherwise
        error('gui:i_ordergrouplevels:badSort', ...
            'Unknown sort order "%s". Use "natural", "count", or "none".', ...
            sortby);
end

names = cLi(ord);
counts = counts(ord);
end
