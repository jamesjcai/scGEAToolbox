function i_invertcolor(src, ~)
    [parentfig] = gui.gui_getfigsce(src);
    axs = findobj(parentfig, 'Type', 'axes');
    if isempty(axs), return; end

    for k = 1:length(axs)
        ax = axs(k);
        cm = colormap(ax);
        if all(cm(1,:) == [.8 .8 .8])
            colormap(ax, [[.8 .8 .8]; flipud(cm(2:end,:))]);
        else
            colormap(ax, flipud(cm));
        end
    end
end
