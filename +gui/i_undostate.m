function i_undostate(src)
%I_UNDOSTATE Label and enable the Undo item to match what is stored.
%
%   gui.i_undostate(app)
%
% Called from EDITMENU's own MenuSelectedFcn, which fires when the menu is
% expanded -- the supported hook for adjusting child items, and the only
% moment the label can be made to name the operation it would reverse.
% Setting Enable and Text here is safe; adding or removing items from a
% menu-expand callback is not, and MathWorks warns it can blank the menu.
%
% This also has to run after IN_ENDISABLEMENU, which turns every
% data-dependent item on once a dataset is loaded and would otherwise leave
% Undo enabled with nothing behind it. Expanding the menu happens later
% than that by definition.
%
% See also GUI.I_SNAPSHOT, GUI.CALLBACK_UNDO.

if nargin < 1 || isempty(src), return; end
if ~isa(src, 'matlab.apps.AppBase'), return; end
if ~isprop(src, 'UndoMenu') || ~pkg.i_isvalid(src.UndoMenu), return; end

m = src.UndoMenu;
s = [];
if pkg.i_isvalid(src.UIFigure)
    s = getappdata(src.UIFigure, 'sceundo');
end

if isempty(s) || ~isstruct(s) || ~isfield(s, 'sce') || isempty(s.sce)
    m.Enable = 'off';
    m.Text = 'Undo';
    return;
end

m.Enable = 'on';
verb = 'Undo';
if isfield(s, 'isredo') && s.isredo
    verb = 'Redo';
end
m.Text = sprintf('%s %s', verb, s.label);

end
