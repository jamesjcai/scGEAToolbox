function callback_ScatterStemPlot(src, ~)

if isa(src,"SingleCellExperiment")
    sce = src;
    FigureHandle = [];
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end

% The plot has no grouping of its own, so narrowing to a subset of cells is
% an opt-in first step; the group chooser is the one the Dotplot, Heatmap
% and Violin plot use - see gui.i_selectgroupsubset.
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
s = sce.s(picked, :);

answer = gui.myQuestdlg(FigureHandle, "Scatter-stem plot for gene expression or cell state variables?","", ...
{'Gene Expression', 'Cell State'}, 'Gene Expression');
switch answer
    case 'Gene Expression'
        [glist] = gui.i_selectngenes(sce, [], FigureHandle);
        if isempty(glist)
            gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
            return;
        end

        [Xt] = gui.i_transformx(sce.X, [], [], FigureHandle);
        if isempty(Xt), return; end

        % answer = gui.myQuestdlg(FigureHandle, 'Plot all in the same figure?','');
        % if strcmp(answer, 'Yes')
        %     fw = gui.myWaitbar(FigureHandle);
        %     gui.i_violinmatrix(full(Xt), sce.g, c, cL, glist, ...
        %             false, '', FigureHandle);
        %
        %     gui.myWaitbar(FigureHandle, fw);
        %
        %     return;
        % end

        n = length(glist);
        thisyv = cell(n,1);
        for k=1:n
            % Subset after normalising, so values match an all-cells plot.
            thisyv{k} = full(Xt(upper(sce.g) == upper(glist(k)), picked));
        end
        ylabelv = glist;

        fw = gui.myWaitbar(FigureHandle);
        gui.sc_uitabgrpfig_feaplot(thisyv, ylabelv, s, FigureHandle, 2);
        gui.myWaitbar(FigureHandle, fw);

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
        else
            gui.myHelpdlg(FigureHandle, 'No valid cell state variables.');
            return;
        end
        for k = 1:length(thisyv)
            thisyv{k} = thisyv{k}(picked);
        end
        gui.sc_uitabgrpfig_feaplot(thisyv, ylabelv, s, FigureHandle, 2);

    otherwise
        return;
end

end
