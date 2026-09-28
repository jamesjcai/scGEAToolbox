function newFigs = i_openinstead(fromFig, openfcn)
%I_OPENINSTEAD Open a window in place of FROMFIG, which comes back after.
%   newFigs = gui.i_openinstead(fromFig, openfcn) calls OPENFCN, and each
%   window it opens takes FROMFIG's place on screen: FROMFIG is hidden, the
%   new window is centred where it was and gets a Back button, and FROMFIG
%   reappears as it was left once the new windows are closed, by Back or
%   otherwise.
%
%   For a button on a table or plot window that opens a window of another
%   kind (a plot from a uifigure table, a table from a plot), where
%   gui.myFigure.drawInto cannot draw into the same window. Without it the
%   new window stacked over the main app and the one it came from.
%
%   Nothing changes when OPENFCN opens no window (it was cancelled, or went
%   to the web browser), or when FROMFIG is not a valid figure.
%
%   See also GUI.MYFIGURE.DRAWINTO.

before = findall(groot, 'Type', 'figure');
openfcn();
newFigs = gobjects(0);
if ~pkg.i_isvalid(fromFig), return; end

after = findall(groot, 'Type', 'figure');
newFigs = after(~ismember(after, before));
% A progress bar left open by OPENFCN is not the window it opened.
isShown = arrayfun(@(f) isvalid(f) && strcmp(f.Visible, 'on') ...
    && ~strcmp(f.Tag, 'TMWWaitbar'), newFigs);
newFigs = newFigs(isShown);
if isempty(newFigs), return; end

for k = 1:numel(newFigs)
    gui.i_movegui2parent(newFigs(k), fromFig);
    in_addback(newFigs(k));
    addlistener(newFigs(k), 'ObjectBeingDestroyed', ...
        @(src, ~) in_restore(src, fromFig, newFigs));
end
fromFig.Visible = 'off';
figure(newFigs(1));
end


function in_addback(fig)
% First on the window's own toolbar, or on a toolbar of its own.
tbs = findall(fig, 'Type', 'uitoolbar', 'Visible', 'on');
custom = tbs(arrayfun(@(t) ~strcmp(t.Tag, 'FigureToolBar'), tbs));
if ~isempty(custom)
    tb = custom(1);
elseif ~isempty(tbs)
    tb = tbs(1);
else
    tb = uitoolbar(fig);
end
gui.i_addbutton2fig(tb, 'off', @(~, ~) close(fig), gui.myFigure.backIcon(), ...
    'Back to the window this was opened from');
% A new tool goes last; move it to the front.
tb.Children = tb.Children([2:end, 1]);
end


function in_restore(closing, fromFig, newFigs)
if ~pkg.i_isvalid(fromFig), return; end
others = newFigs(isvalid(newFigs) & newFigs ~= closing);
others = others(arrayfun(@(f) strcmp(f.BeingDeleted, 'off'), others));
if isempty(others)
    fromFig.Visible = 'on';
    figure(fromFig);
end
end
