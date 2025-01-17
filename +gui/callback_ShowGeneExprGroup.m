function callback_ShowGeneExprGroup(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    % [axx, bxx] = view(findall(FigureHandle,'type','axes'));    

    allowunique = false;
    [thisc] = gui.i_select1class(sce, allowunique);
    if isempty(thisc), return; end

    answer = questdlg("Select a dependent variable from gene expression or cell state?","", ...
        'Gene Expression', 'Cell State','Gene Expression');

    switch answer
        case 'Gene Expression'
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist), return; end
            % answer = questdlg("Select the type of expression values","",...
            %     "Raw UMI Counts","Library Size-Normalized",)            

            gui.i_feaplotarray(sce, glist, thisc, false, FigureHandle);
        case 'Cell State'
            a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
            s = sce.s;
            [thisyv, ylabelv] = gui.i_selectnstates(sce, true);
            y = thisyv{1};
            ylabelv = string(ylabelv);

            [c, cL, noanswer] = gui.i_reordergroups(thisc);
            if noanswer, return; end

            hx = gui.myFigure;
            for ky = 1:length(cL)
                cellidx = c==ky;
                nexttile;
                ydata = y(cellidx);
                if size(s,2)>2
                    scatter3(s(cellidx, 1), s(cellidx, 2), s(cellidx, 3), 5, ydata, 'filled');
                else
                    scatter(s(cellidx, 1), s(cellidx, 2), 5, ydata, 'filled');
                end
                gui.i_setautumncolor(ydata, a, true, any(ydata==0));
                z = ydata;
                clim([min(z) max(z)]);  % Adjust color axis to data range
                title(cL{ky});
            end
            sgtitle(strrep(ylabelv,'_','\_'));
            hx.show(FigureHandle);
    end
            
end