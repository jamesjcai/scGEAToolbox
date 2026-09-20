function varargout = scgeatool(varargin)
% SCGEATOOL - Launch App Designer GUI
%
% Use:
%   scgeatool            % launches App Designer GUI
%
% Every menu item and the whole menu arrangement are baked into
% scgeatoolApp.mlapp, so this launches the app and nothing else.
%
% It used to call a series of functions that shaped the menus at launch,
% because a binary .mlapp was thought to be uneditable. Seven of them added
% a menu item; those are deleted, the items being in the .mlapp and the
% functions unreachable from anything but their own tests. What is kept is
% the half that reorganises and relabels:
%
%   shape   - i_mergeharmonymenu, i_regroupexternalmenu, i_groupeditmenu,
%             i_groupanalyzemenu, i_groupannotatemenu (all on i_groupmenu)
%   labels  - i_addmenuicons, i_addmenutooltips
%
% Their tables are the written record of why each menu is arranged and
% labelled as it is, which is worth more than the code around them, and
% each is idempotent, so calling one against the app is a no-op.
% TESTS/MENUSHAPERTEST is what holds that true. It had already stopped
% being true once, and quietly: nothing about a duplicated submenu or an
% item sinking to the bottom of its menu raises anything.
%
% Calling them is no longer needed, and doing it at launch would also miss
% the eight places that construct SCGEATOOLAPP directly -- which is the
% reason the arrangement is baked in rather than injected here.

if ~gui.i_installed('stats'), return; end
app = scgeatoolApp(varargin{:});
if nargout > 0
    varargout{1} = app;
end
end
