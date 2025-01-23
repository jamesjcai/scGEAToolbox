function i_pickcolor(~, ~, revcolor, hFig)

if nargin < 4, hFig = gcf; end
if nargin < 3, revcolor = false; end
list = {
    'parula', 'jet', 'hsv', 'hot', 'cool', ...
    'spring', 'summer', 'autumn', 'winter', ...
    'gray', 'bone', 'copper', 'pink', 'lines', ...
    'colorcube', 'prism', 'flag', 'white', 'turbo'};
[indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
    'PromptString', 'Select a colormap:');
if tf == 1
    a = colormap(list{indx});
    if revcolor
        colormap(hFig, flipud(a));
    else
        colormap(hFig, a);
    end
end
end