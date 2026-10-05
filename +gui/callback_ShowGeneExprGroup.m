function callback_ShowGeneExprGroup(src, ~)


  [FigureHandle, sce] = gui.gui_getfigsce(src);

    allowunique = false;
    % Several grouping variables may be picked; they cross into one composite
    % label per cell ("Macrophages | IL"). Downstream treats thisc as a
    % per-cell label vector, so the composite needs no special handling.
    [thisc, clabel] = gui.i_selectnclass(sce, allowunique,[],[],FigureHandle);
    if isempty(thisc), return; end
    thisc = string(thisc);

    % Same group chooser as gui.callback_Dotplot - see gui.i_selectgroupsubset.
    [picked, levels] = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
    if isempty(picked), return; end
    thisc = thisc(picked);

    answer = gui.myQuestdlg(FigureHandle, ...
        "Select a dependent variable from gene expression or cell state?","", ...
        {'Gene Expression', 'Cell State'},'Gene Expression');

    switch answer
        case 'Gene Expression'
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist), return; end
            % answer = gui.myQuestdlg(FigureHandle, "Select the type of expression values","",...
            %     "Raw UMI Counts","Library Size-Normalized",)

            gui.i_feaplotarray(sce, glist, thisc, false, FigureHandle, levels, picked);
        case 'Cell State'
            a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
            s = sce.s(picked, :);
            [thisyv, ylabelv] = gui.i_selectnstates(sce, true, [1], FigureHandle);
            if isempty(thisyv) || isempty(ylabelv), return; end

            y = thisyv{1}(picked);
            ylabelv = string(ylabelv);

            [c, cL, noanswer] = gui.i_reordergroups(thisc, levels, FigureHandle);
            if noanswer, return; end

            % One colour scale for every panel, as the gene path has, so the
            % groups can be compared. Each panel used to take its own range,
            % and a group whose values were all equal gave CLIM two equal
            % limits, which it rejects.
            yall = y(c > 0);
            ylim2 = [min(yall), max(yall)];
            hx = gui.myFigure(FigureHandle);
            for ky = 1:length(cL)
                cellidx = c==ky;
                nexttile;
                ydata = y(cellidx);
                if size(s,2)>2
                    scatter3(s(cellidx, 1), s(cellidx, 2), s(cellidx, 3), 5, ydata, 'filled');
                else
                    scatter(s(cellidx, 1), s(cellidx, 2), 5, ydata, 'filled');
                end
                gui.i_setautumncolor(ydata, a, true, any(ydata==0), [], FigureHandle);
                if ylim2(2) > ylim2(1)
                    clim(ylim2);
                end
                title(cL{ky});
            end
            sgtitle(strrep(ylabelv,'_','\_'));
            hx.show(FigureHandle);
    end

end
