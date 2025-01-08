function i_twogenestemscatter(sce, g1, g2)

if ismember(upper(g1), upper(sce.g)) && ismember(upper(g2), upper(sce.g))

    hx=gui.myFigure;
    hFig=hx.FigureHandle;



    subplot(1, 2, 1);
    sc_scattermarker(sce.X, upper(sce.g), sce.s, ...
        upper(g1), 1, [], false);
    box on
    subplot(1, 2, 2);
    sc_scattermarker(sce.X, upper(sce.g), sce.s, ...
        upper(g2), 1, [], false);
    box on
    hFig.Position(3) = hFig.Position(3) * 1.55;
    hx.show;
else
    disp('aaa');
end