function s = i_effectiveundo(fig, livesce)
%I_EFFECTIVEUNDO The snapshot Undo should actually restore.
%
%   s = gui.i_effectiveundo(app.UIFigure, app.sce)
%
% GUI.I_SNAPSHOT runs at the head of every destructive handler, before
% its dialogs, so a handler the user then cancels still overwrites the
% snapshot -- with a copy of data it never changed. Undo would then
% restore that unchanged copy and the operation before it could no longer
% be undone: delete cells, open Filter and cancel, and the deleted cells
% were gone for good.
%
% I_SNAPSHOT therefore keeps the snapshot it replaces one level down, in
% PREV. When the live data still equals the newest snapshot, the operation
% after it changed nothing, and this returns PREV instead. ISEQUAL
% compares property values, display state included, and costs about a
% millisecond on an unchanged copy because COPY shares the arrays.
%
% Returns [] when there is nothing to undo.
%
% See also GUI.I_SNAPSHOT, GUI.CALLBACK_UNDO, GUI.I_UNDOSTATE.

s = [];
if isempty(fig) || ~pkg.i_isvalid(fig), return; end
s = getappdata(fig, 'sceundo');
if isempty(s) || ~isstruct(s) || ~isfield(s, 'sce') || isempty(s.sce)
    s = [];
    return;
end
if ~s.isredo && ~isempty(livesce) && isequal(livesce, s.sce)
    % The newest operation changed nothing (it was cancelled).
    if isfield(s, 'prev') && ~isempty(s.prev)
        s = s.prev;
    else
        s = [];
    end
end
end
