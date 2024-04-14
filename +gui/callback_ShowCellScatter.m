function callback_ShowCellScatter(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));

    [thisx, xlabelv] = gui.i_select1state(sce, false, false, false);
    if isempty(thisx), return; end


    answer = questdlg("Violinplot for gene expression or cell state variables?","", ...
        'Gene Expression', 'Cell State','Gene Expression');

    switch answer
        case 'Gene Expression'
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist)
                helpdlg('No gene selected.', '');
                return;
            end
            gui.sc_uitabgrpfig_scatter(sce, glist, thisx, xlabelv, FigureHandle);
        case 'Cell State'
            [thisyv, ylabelv] = gui.i_selectnstates(sce, true);
            a = false(length(thisyv), 1);
            for k = 1:length(thisyv)
                a(k) = isnumeric(thisyv{k});
            end
            if any(a)
                if ~all(a)
                    thisyv = thisyv(a);
                    ylabelv = ylabelv(a);
                    waitfor(helpdlg('Only continuous variables of cell state will be shown.',''));
                end
                gui.i_scattertabs(thisyv, ylabelv, thisx, xlabelv, FigureHandle);
            else
                waitfor(helpdlg('No valid cell state variables. Violinplot cannot be shown.',''));
            end            
    end
end
    
