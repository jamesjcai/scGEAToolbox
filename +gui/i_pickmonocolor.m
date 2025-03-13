function i_pickmonocolor(~, ~, revcolor, parentfig)

if nargin < 4, parentfig = gcf; end
if nargin < 3, revcolor = true; end
list = {'spring', ...
    'summer', 'autumn', ...
    'gray', 'bone', 'pink', 'copper'};
    if gui.i_isuifig(parentfig)
        [indx, tf] = gui.myListdlg(parentfig, list, 'Select a colormap:');
    else
        [indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
            'PromptString', 'Select a colormap:');
    end
if tf == 1
    a = colormap(list{indx});
    if revcolor
        colormap(parentfig, flipud(a));
    else
        colormap(parentfig, a);
    end
end
end