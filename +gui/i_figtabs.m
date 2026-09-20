function tabs = i_figtabs(fig)
%I_FIGTABS The tabs of a figure's tab group, in the order they are shown.
%
%   tabs = gui.i_figtabs(fig)
%
% Returns the tabs of FIG's tab group, left to right, or an empty GOBJECTS
% array when FIG has no tab group - or more than one.
%
% More than one is deliberately not answered. Two groups side by side have no
% single list of tabs: what is on screen is one selection from each, so
% "every tab" would mean a slide per combination. A nested group has the same
% problem from the other end, since its parent tab is in every shot. A caller
% that gets nothing back should treat the figure as the single view it
% appears to be, which is what it did before tabs were considered at all.
%
% see also: gui.i_export2pptx, gui.myFigure

tabs = gobjects(0);
if nargin < 1 || isempty(fig) || ~pkg.i_isvalid(fig), return; end

tg = findall(fig, 'Type', 'uitabgroup');
if ~isscalar(tg), return; end

% Children, not FINDALL: a tab group reports its tabs in the order it draws
% them, so slide 1 is the leftmost tab rather than whichever was made last.
tabs = tg.Children;
tabs = tabs(arrayfun(@(h) isa(h, 'matlab.ui.container.Tab'), tabs));
end
