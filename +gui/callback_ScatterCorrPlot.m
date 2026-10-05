function callback_ScatterCorrPlot(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, 'Select an independent variable. Continue?','');
if ~strcmp(answer,'Yes'), return; end

% The plot has no grouping of its own, so narrowing to a subset of cells is
% an opt-in step; the group chooser is the one the Dotplot, Heatmap and
% violin plots use - see gui.i_selectgroupsubset.
picked = true(size(sce.X, 2), 1);
answer = gui.myQuestdlg(FigureHandle, ...
    "Plot all cells, or only cells in selected groups?", "", ...
    {'All Cells', 'Selected Groups', 'Cancel'}, 'All Cells');
switch answer
    case 'All Cells'
        % Keep every cell.
    case 'Selected Groups'
        [thisc, clabel] = gui.i_selectnclass(sce, false, [], [], FigureHandle);
        if isempty(thisc), return; end
        picked = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
        if isempty(picked), return; end
    otherwise
        return;
end

[thisx, xlabelv] = gui.i_select1state(sce, false, false, false, true, FigureHandle);
if isempty(thisx), return; end
if ~isnumeric(thisx)
    gui.myWarndlg(FigureHandle, 'This function works with continuous variables only.');
    return;
end
thisx = thisx(picked);

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
            % Subset after normalising, so values match an all-cells plot.
            y{k} = full(Xt(upper(sce.g) == upper(glist(k)), picked));
        end
        gui.i_scattertabs(y, glist, thisx, xlabelv, FigureHandle);


       % i_plot_pseudotimeseries(X, genelist, t, genes)
       % %Plot pseudotime series

    case 'Cell State'
        [thisyv, ylabelv] = gui.i_selectnstates(sce, true, [1], FigureHandle);
        % Cancelled, or I_SELECTNSTATES has already said there is nothing
        % to pick: return quietly rather than report "No valid cell state
        % variables", which is what an empty pick used to fall through to.
        if isempty(thisyv), return; end
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
            for k = 1:length(thisyv)
                thisyv{k} = thisyv{k}(picked);
            end
            gui.i_scattertabs(thisyv, ylabelv, thisx, xlabelv, FigureHandle);
        else
            gui.myHelpdlg(FigureHandle, 'No valid cell state variables.');
        end
end

end
