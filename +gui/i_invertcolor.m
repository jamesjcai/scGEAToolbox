function i_invertcolor(~, ~)
    cm = colormap();
    if all(cm(1,:) == [.8 .8 .8])
        colormap([[.8 .8 .8]; flipud(cm(2:end,:))]);
    else
        colormap(flipud(cm));
    end
end
