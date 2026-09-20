function [ok] = i_setbrushdata(h, mask)
%I_SETBRUSHDATA Make the plot highlight the given cells.
%
%   ok = gui.i_setbrushdata(h, mask)
%
%   H is a scatter handle, MASK a logical selection with one entry per point.
%   Returns false, and changes nothing, when H is not a usable scatter or when
%   MASK is not the length of its data - a plot of something other than these
%   cells, where writing the mask would highlight arbitrary points.
%
%   BRUSHDATA is a uint8 row the same shape as the plotted data, and it is
%   what the data brush and every handler that reads a selection go by, so
%   writing it is how a selection worked out in code becomes the selection the
%   user can see.
%
%   See also gui.i_expandbrushed, gui.callback_ExpandBrushedCells.

ok = false;
if isempty(h) || ~pkg.i_isvalid(h) || ~isprop(h, 'BrushData'), return; end
if isempty(h.BrushData) || numel(h.BrushData) ~= numel(mask), return; end

h.BrushData = uint8(reshape(logical(mask), size(h.BrushData)));
ok = true;
end
