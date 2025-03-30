function i_resizewin(~, ~, hFig)
    w = gui.i_inputnumk(450, 10, 2000, 'Window width', hFig);
    if isempty(w), return; end
    h = gui.i_inputnumk(420, 10, 2000, 'Window height', hFig);
    if isempty(h), return; end
    hFig.Position = [hFig.Position(1) hFig.Position(2) w h];
end