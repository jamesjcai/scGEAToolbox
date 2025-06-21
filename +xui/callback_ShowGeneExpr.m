function callback_ShowGeneExpr(app, ~)
    [glist] = gui.i_selectngenes(app.sce, [], app.UIFigure);
    if isempty(glist), return; end
    [Xt] = gui.i_transformx(app.sce.X, [], [], app.UIFigure);
    if isempty(Xt), return; end
    n = length(glist);
    for k = 1:n
        figure;
        sc_scattermarker(Xt, app.sce.g, app.sce.s, glist(k), 2, 5, false);
    end
end