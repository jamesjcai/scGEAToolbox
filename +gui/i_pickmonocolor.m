function i_pickmonocolor(~, ~, revcolor, hFig)

if nargin < 4, hFig = gcf; end
if nargin < 3, revcolor = true; end
list = {'spring', ...
    'summer', 'autumn', ...
    'gray', 'bone', 'pink', 'copper'};
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