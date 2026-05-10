function i_movegui2parent(hFig, parentfig)
try
    dlgSize = hFig.Position(3:4);
    pos = gui.i_centerdlgpos(parentfig, dlgSize);
    hFig.Position = pos;
catch
end
end
