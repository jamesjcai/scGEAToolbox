function [newc, p, neworder] = i_grouporderperm(c, cL, ordnow, sortby)
%I_GROUPORDERPERM Reorder the group blocks of an already-grouped plot.
%
%   [newc, p, neworder] = gui.i_grouporderperm(c, cL, ordnow, sortby)
%
%   Inputs:
%     c       - per-cell group index, 1..K, already sorted so that the
%               cells of a group sit together.
%     cL      - the K group labels, in the order they are drawn in.
%     ordnow  - which original group each drawn position holds, so that
%               "none" has an order to go back to. (1:K).' before any sort.
%     sortby  - "count" (most cells first), "natural" (alphabetic, with
%               embedded numbers read as numbers), or "none" (the order
%               the groups arrived in).
%
%   Outputs:
%     newc     - c renumbered into the new order, and re-sorted.
%     p        - the permutation of the drawn positions. Apply it to the
%                labels and to ORDNOW: cL(p), ordnow(p).
%     neworder - the new column order. Apply it to the data: Y(:, neworder).
%
%   The permutation is worked out against the order on screen and applied
%   to it, rather than rebuilt from the original each time, which is what
%   lets a re-sort compose with a rename.
%
%   See also GUI.I_ORDERGROUPLEVELS, GUI.I_ASKGROUPORDER, GUI.I_HEATMAP.

switch string(sortby)
    case "count"
        % SORT is stable, so equal counts keep the order they are already
        % in rather than jumping about on every click.
        [~, p] = sort(accumarray(c(:), 1), 'descend');
    case "natural"
        [~, p] = natsort(string(cL));
    case "none"
        [~, p] = sort(ordnow);
    otherwise
        error('gui:i_grouporderperm:badSort', ...
            'Unknown sort order "%s". Use "count", "natural", or "none".', ...
            string(sortby));
end
p = p(:);

% Renumber the groups into the new order, then sort the cells by the new
% number. SORT is stable, so cells keep their relative order inside a group.
key = zeros(size(c));
for kg = 1:numel(p)
    key(c == p(kg)) = kg;
end
[newc, neworder] = sort(key);
end
