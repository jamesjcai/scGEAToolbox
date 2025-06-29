function callback_ScatterCorrPlot(src, ~)


        [FigureHandle, sce] = gui.gui_getfigsce(src);

    %[axx, bxx] = view(findall(FigureHandle,'type','axes'));

    answer = gui.myQuestdlg(FigureHandle, 'Select an independent variable. Continue?','');
    if ~strcmp(answer,'Yes'), return; end

    [thisx, xlabelv] = gui.i_select1state(sce, false, false, false, true, FigureHandle);
    if isempty(thisx), return; end
    if ~isnumeric(thisx)
        gui.myWarndlg(FigureHandle, 'This function works with continuous varibles only.');
        return;
    end

    answer = gui.myQuestdlg(FigureHandle, "Select a dependent variable from gene expression or cell state?","", ...
        {'Gene Expression', 'Cell State'},'Gene Expression');

    switch answer
        case 'Gene Expression'
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist)
                gui.myHelpdlg(FigureHandle, 'No gene selected.');
                return;
            end
            [Xt] = gui.i_transformx(sce.X, [], [], FigureHandle);
            if isempty(Xt), return; end
            n = length(glist);
            y=cell(n,1);
            for k=1:n
                y{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
            end
            gui.i_scattertabs(y, glist, thisx, xlabelv, FigureHandle);


       
           % i_plot_pseudotimeseries(X, genelist, t, genes)
           % %Plot pseudotime series

        case 'Cell State'
            [thisyv, ylabelv] = gui.i_selectnstates(sce, true, [], FigureHandle);
            a = false(length(thisyv), 1);
            for k = 1:length(thisyv)
                a(k) = isnumeric(thisyv{k});
            end
            if any(a)
                if ~all(a)
                    thisyv = thisyv(a);
                    ylabelv = ylabelv(a);
                    gui.myHelpdlg(FigureHandle, 'Only continuous variables of cell state will be shown.');
                end
                gui.i_scattertabs(thisyv, ylabelv, thisx, xlabelv, FigureHandle);
            else
                gui.myHelpdlg(FigureHandle, 'No valid cell state variables.');
            end            
    end
    
end    
