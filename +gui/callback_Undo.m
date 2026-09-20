function [done] = callback_Undo(src, ~)
%CALLBACK_UNDO Put the dataset back the way GUI.I_SNAPSHOT last saw it.
%
%   gui.callback_Undo(app)
%
% Swaps rather than restores: the state being replaced becomes the new
% snapshot, so the same menu item undoes and then redoes. That costs no
% more than restoring would -- the object being swapped out already exists
% -- and it means a mis-aimed Undo is itself undoable.
%
% See also GUI.I_SNAPSHOT, GUI.I_UNDOSTATE.

done = false;
if nargin < 1 || isempty(src), return; end
if ~isa(src, 'matlab.apps.AppBase'), return; end
if ~pkg.i_isvalid(src.UIFigure), return; end

s = getappdata(src.UIFigure, 'sceundo');
if isempty(s) || ~isstruct(s) || ~isfield(s, 'sce') || isempty(s.sce)
    gui.myHelpdlg(src.UIFigure, 'Nothing to undo.');
    return;
end

% The outgoing state becomes the way back. COPY is not needed: the live
% object is about to be replaced rather than mutated, so handing the old
% handle to the snapshot is safe and saves a duplicate.
previous = src.sce;
src.sce = s.sce;

s.sce = previous;
s.isredo = ~s.isredo;
setappdata(src.UIFigure, 'sceundo', s);

[src.c, src.cL] = findgroups(string(src.sce.c));
src.in_RefreshAll(true, false);
done = true;

end
