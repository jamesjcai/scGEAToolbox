function callback_PickMarkerType(src, ~)
%CALLBACK_PICKMARKERTYPE  Choose the scatter plot marker symbol.
%
%   gui.callback_PickMarkerType(app)
%
%   Lists the marker symbols and applies the chosen one to the current
%   scatter plot. The toolbar button sharing this feature steps to the next
%   marker without asking; this is the entry that lets one be picked, the
%   same split GUI.CALLBACK_PICKCOLORMAP makes between the colormap button
%   and the colormap menu item.
%
%   A point marker needs a far larger SizeData than a shaped one to stay
%   visible, so the size is adjusted when crossing between the two. Picking
%   one shape over another leaves the size alone.
%
% See also GUI.CALLBACK_PICKCOLORMAP, GUI.I_GSCATTER3.

if isa(src, 'matlab.apps.AppBase')
    parentfig = src.UIFigure;
    h = src.h;
else
    % A figure is handled directly, as GUI.CALLBACK_SAVEX does:
    % GUI_GETFIGSCE reads src.Parent.Parent, which is the root for a
    % figure and leaves it with nothing to return.
    if isa(src, 'matlab.ui.Figure')
        parentfig = src;
    else
        parentfig = gui.gui_getfigsce(src);
    end
    h = findobj(parentfig, 'Type', 'scatter');
end

if isempty(h) || ~all(pkg.i_isvalid(h))
    gui.myWarndlg(parentfig, 'No scatter plot available.');
    return;
end

% MATLAB's scatter markers, least to most ornate. '.' is first because it
% is the toolbox default and the one most plots start from.
% A column, like NAMES below: with one a row and the other a column, the
% two comparisons further down expand into a 13-by-13 matrix and FIND
% returns an index off the end of both.
symbols = [".", "o", "+", "*", "x", "s", "d", "^", "v", ">", "<", "p", "h"]';
names = ["point  ."
         "circle  o"
         "plus  +"
         "asterisk  *"
         "cross  x"
         "square  s"
         "diamond  d"
         "triangle up  ^"
         "triangle down  v"
         "triangle right  >"
         "triangle left  <"
         "pentagram  p"
         "hexagram  h"];

% MATLAB hands some markers back under their long name - 's' reads as
% "square" - so the current one is matched against both spellings and
% normalized to its symbol before anything below compares it.
current = string(h(1).Marker);
descriptions = strtrim(extractBefore(names + "  ", "  "));
k = find(symbols == current | descriptions == current, 1);
prefer = [];
if ~isempty(k)
    prefer = char(names(k));
    current = symbols(k);
end

[indx, tf] = gui.myListdlg(parentfig, cellstr(names), ...
    'Pick a marker type:', prefer, false);
if tf ~= 1 || isempty(indx), return; end

chosen = symbols(indx);
set(h, 'Marker', char(chosen));

% Only scatter objects carry SizeData; GSCATTER returns lines, which size
% themselves through MarkerSize instead.
if ~isprop(h(1), 'SizeData')
    return;
end
if chosen == "." && current ~= "."
    set(h, 'SizeData', 50);
elseif chosen ~= "." && current == "."
    set(h, 'SizeData', 10);
end

% Record it on the SCE so the choice survives a save and reload.
gui.i_storedisplay(src);
end
