function i_pickmonocolor(~, ~, revcolor)
if nargin < 3, revcolor = true; end
list = {'spring', ...
    'summer', 'autumn', ...
    'gray', 'bone', 'pink', 'copper'};
[indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
    'PromptString', 'Select a colormap:');
if tf == 1
    a = colormap(list{indx});
    if revcolor
        colormap(flipud(a));
    else
        colormap(a);
    end
end
end