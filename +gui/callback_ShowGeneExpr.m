function callback_ShowGeneExpr(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end

    % answer = questdlg("Select the type of expression values","",...
    %     "Raw UMI Counts","Library Size-Normalized",)

    [Xt] = gui.i_transformx(sce.X);
    if isempty(Xt), return; end

    answer = questdlg("Plot on the same figure?","", "Yes, same figure", ...
        "No, different figures", "Cancel", "Yes, same figure");
    if strcmp(answer, "Cancel"), return; end

    fw=gui.gui_waitbar;
    n = length(glist);

    switch answer
        case "No, different figures"
            y = cell(n,1);
            for k = 1:n
                y{k} = Xt(sce.g == glist(k), :);
            end
            gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
        case "Yes, same figure"
            hFig = figure(Visible="off");
            maxy = 0;
            for k = 1:n
                nexttile
                sc_scattermarker(Xt, sce.g, sce.s, glist(k), 2, 5, false);
                colorbar;
                maxy = max([maxy, max(Xt(sce.g == glist(k)))]);
            end
            gui.i_movegui2parent(hFig, FigureHandle);
            % tb = uitoolbar(hFig);
            % uipushtool(tb, 'Separator', 'off');
            % pkg.i_addbutton2fig(tb, 'off', {@in_rescale, maxy}, 'networkcomp.gif', 'Normalize scales...');
            hFig.Visible = "on";
    end
    gui.gui_waitbar(fw);

    
    % function in_rescale(~, ~, maxy)
    %        clim([0 maxy]);
    % end

end


