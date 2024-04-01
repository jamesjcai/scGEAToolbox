function i_movegui2parent(hFig, parentfig)
    
    try
        if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
            [px_new] = gui.i_getchildpos(parentfig, hFig);
            if ~isempty(px_new)
                movegui(hFig, px_new);
            else
                movegui(hFig, 'center');
            end
        else
            movegui(hFig, 'center');
        end
    catch
        movegui(hFig, 'center');
    end
end