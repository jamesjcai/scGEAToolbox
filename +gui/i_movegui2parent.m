function i_movegui2parent(hFig, parentfig)
try
    dlgSize = hFig.Position(3:4);
    pos = gui.i_centerdlgpos(parentfig, dlgSize);
    hFig.Position = pos;
catch
    % dialog stays at its default position if centering fails (no parent etc.)
end
end
