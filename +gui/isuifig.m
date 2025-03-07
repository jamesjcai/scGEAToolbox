function isUIFig = isuifig(fig)
    if isprop(fig, 'Scrollable')
        isUIFig = true;
    else
        isUIFig = false;
    end
end


