function i_pickcolormap(~, ~, ~, setzerocolor)
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
[indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
    'PromptString', 'Select a colormap:');
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