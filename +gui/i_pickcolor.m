function i_pickcolor(src, ~, revcolor, parentfig)

% The toolbar button's own figure is the target. gcf would miss a uifigure
% (HandleVisibility is off) and land on whichever figure was current.
if nargin < 4 || isempty(parentfig)
    parentfig = ancestor(src, 'figure');
    if isempty(parentfig), parentfig = gcf; end
end
if nargin < 3, revcolor = false; end
cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
list = {
'parula', 'jet', 'hsv', 'hot', 'cool', ...
'spring', 'summer', 'autumn', 'winter', ...
'gray', 'bone', 'copper', 'pink', 'lines', ...
'colorcube', 'prism', 'flag', 'white', 'turbo'};
if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, list, 'Select a colormap:', [], false);
else
    [indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
        'PromptString', 'Select a colormap:');
end
if tf == 1
    % Build the map directly; colormap(name) would apply it to gcf.
    a = feval(list{indx}, size(colormap(parentfig), 1));
    if revcolor
        a = flipud(a);
    end
    % A figure colormap resets every axes in it, those with their own map
    % (Invert Colors sets one per axes) included -- measured on classic
    % and ui figures -- so this one call covers all the panels.
    colormap(parentfig, a);
end
end
