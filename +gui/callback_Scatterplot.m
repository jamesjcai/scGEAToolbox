function callback_Scatterplot(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    %[axx, bxx] = view(findall(FigureHandle,'type','axes'));

    [thisx, xlabelv] = gui.i_select1state(sce, false, false, false, true);
    if isempty(thisx), return; end
    if ~isnumeric(thisx)
        warndlg('This function works with continuous varibles only.','');
        return;
    end



    answer = questdlg("Violinplot for gene expression or cell state variables?","", ...
        'Gene Expression', 'Cell State','Gene Expression');

    switch answer
        case 'Gene Expression'
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist)
                helpdlg('No gene selected.', '');
                return;
            end
            [Xt] = gui.i_transformx(sce.X);
            n = length(glist);
            y=cell(n,1);
            for k=1:n
                y{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
            end
            gui.i_scattertabs(y, glist, thisx, xlabelv, FigureHandle);


            % gui.sc_uitabgrpfig_scatter(sce, glist, thisx, xlabelv, FigureHandle);
        
           % i_plot_pseudotimeseries(X, genelist, t, genes)
           % %Plot pseudotime series

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
                waitfor(helpdlg('No valid cell state variables.',''));
            end            
    end
end
    
