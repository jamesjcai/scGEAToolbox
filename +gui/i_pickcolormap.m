function i_pickcolormap(src, ~, ~, setzerocolor, parentfig)

% The clicked button's figure is the target: gcf never returns a uifigure.
if nargin < 5 || isempty(parentfig), parentfig = ancestor(src, 'figure'); end
if nargin < 4, setzerocolor = false; end
if nargin < 3, c = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
% A heatmap chart keeps its own Colormap, so the figure's map would miss it.
target = parentfig.CurrentAxes;
if isempty(target)
    target = findobj(parentfig, 'Type', 'axes', '-or', 'Type', 'heatmap');
    if isempty(target), target = parentfig; else, target = target(1); end
end
a = colormap(target);
if all(a(1, :) == [0.8 0.8 0.8])
    setzerocolor = true;
end
if setzerocolor
    list = {'parula', 'turbo', 'hsv', 'hot', 'cool', 'spring', ...
        'summer', 'autumn (default)', ...
        'winter', 'jet'};
else
    list = {'parula', 'turbo', 'hsv', 'hot', 'cool', 'spring', ...
        'summer', 'autumn', ...
        'winter', 'jet', 'gray', 'bone', 'pink', 'copper', 'lines'};
end
if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, list, 'Select a colormap:', [], false);
else
    [indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
        'PromptString', 'Select a colormap:');
end
if tf == 1
    a = list{indx};
    if strcmp(a, 'autumn (default)')
        a = 'autumn';
    end
    if setzerocolor
        gui.i_setautumncolor(1:5, a, [], [], target, parentfig);
    else
        colormap(target, a);
    end
end
end
