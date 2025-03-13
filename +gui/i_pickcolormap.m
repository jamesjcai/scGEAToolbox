function i_pickcolormap(~, ~, ~, setzerocolor, parentfig)

if nargin < 5, parentfig = []; end
if nargin < 4, setzerocolor = false; end
if nargin < 3, c = []; end
a = colormap;
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
    [indx, tf] = gui.myListdlg(parentfig, list, 'Select a colormap:');
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
        gui.i_setautumncolor(1:5, a);
    else
        colormap(a);
    end
    % setpref('scgeatoolbox','prefcolormapname',a);
end
end