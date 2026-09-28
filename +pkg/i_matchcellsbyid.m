function [order, status, msg] = i_matchcellsbyid(targetid, sourceid)
%I_MATCHCELLSBYID Line incoming per-cell values up with this dataset's cells.
%
%   [order, status, msg] = pkg.i_matchcellsbyid(targetid, sourceid)
%
%   targetid  cell IDs of the cells being written to, one per cell
%   sourceid  cell IDs of the incoming object, one per incoming value
%
%   Returns ORDER, indexing the incoming values so that value(ORDER(j))
%   belongs to targetid(j). It is (1:n)' when the two are already in step,
%   and [] when they cannot be lined up at all.
%
%   STATUS says what the IDs showed:
%     "positional"  one side carries no IDs, so order is the only evidence
%     "inorder"     the IDs agree cell for cell
%     "reordered"   the same cells in a different order; ORDER puts them in step
%     "mismatch"    not the same set of cells, or not the same number of them
%     "duplicates"  an ID repeats, so no mapping can be trusted
%
%   MSG is a sentence a caller can show the user, phrased the same way for
%   every status so it can be reported either way.
%
%   Why this exists
%   ---------------
%   A matching cell count is not a matching set of cells. Two SCEs of the
%   same size are routinely two different subsets of the same experiment, and
%   copying one's per-cell values onto the other by column order then
%   mislabels every cell without anything going wrong out loud. Cell IDs are
%   the only evidence that they are the same cells.
%
%   A "mismatch" is not by itself a reason to refuse: cell IDs pick up batch
%   suffixes when SCEs are merged and differ between readers, so IDs that do
%   not line up are sometimes the same cells spelled differently. This
%   function reports; the caller decides.
%
%   Duplicates are different. ISMEMBER resolves a repeated ID by silently
%   taking the first match, which copies one cell's value onto another, so
%   they are reported as their own status rather than folded into a mapping
%   that cannot be right.
%
%   See also GUI.CALLBACK_MERGECELLSUBTYPES, PKG.I_ALIGNMODALITY.

targetid = string(targetid(:));
sourceid = string(sourceid(:));
n = numel(targetid);
order = (1:n)';

% An SCE built from a bare count matrix has no cell IDs, and a list of empty
% strings is that same absence spelled differently.
if in_isblank(targetid) || in_isblank(sourceid)
    status = "positional";
    if in_isblank(targetid)
        msg = "This dataset carries no cell IDs, so the cells can only be matched by position.";
    else
        msg = "The incoming data carries no cell IDs, so the cells can only be matched by position.";
    end
    return;
end

if numel(sourceid) ~= n
    status = "mismatch";
    order = [];
    msg = string(sprintf('The incoming data has %s for %s.', ...
        pkg.i_plural(numel(sourceid), 'cell ID'), ...
        pkg.i_plural(n, 'cell')));
    return;
end

if isequal(targetid, sourceid)
    status = "inorder";
    msg = string(sprintf('All %s matched, in the same order.', ...
        pkg.i_plural(n, 'cell ID')));
    return;
end

dup = in_firstduplicate(targetid);
side = "this dataset";
if ismissing(dup)
    dup = in_firstduplicate(sourceid);
    side = "the incoming data";
end
if ~ismissing(dup)
    status = "duplicates";
    order = [];
    msg = string(sprintf(['The cell ID ''%s'' appears more than once in %s, ', ...
        'so the cells cannot be matched by ID.'], dup, side));
    return;
end

[tf, loc] = ismember(targetid, sourceid);
if all(tf)
    status = "reordered";
    order = loc;
    msg = string(sprintf(['All %s matched, in a different order; the ', ...
        'incoming values were lined up with them.'], ...
        pkg.i_plural(n, 'cell ID')));
else
    status = "mismatch";
    order = [];
    msg = string(sprintf(['%d of the %s are not in the incoming data (for ', ...
        'example ''%s'').'], sum(~tf), pkg.i_plural(n, 'cell ID'), ...
        targetid(find(~tf, 1))));
end
end

function tf = in_isblank(id)
tf = isempty(id) || all(ismissing(id) | strlength(id) == 0);
end

function id = in_firstduplicate(ids)
% The first repeated value, or <missing> when every value is distinct. Named
% so the caller can say which ID is the problem: "an ID repeats" sends the
% reader looking through the whole list.

id = string(missing);
[sorted, ix] = sort(ids);
rep = find([false; sorted(2:end) == sorted(1:end-1)], 1);
if ~isempty(rep)
    id = ids(min(ix(rep), ix(rep-1)));
end
end
