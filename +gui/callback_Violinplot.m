function callback_Violinplot(src, ~)

if isa(src,"SingleCellExperiment")
    sce = src;
    FigureHandle = [];
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end

allowunique = false;
% Several grouping variables may be picked; they are crossed into one
% composite label per cell ("Macrophages | IL"), matching what the cell
% score comparison does. Everything below treats thisc as a per-cell label
% vector, so the composite needs no special handling past this point.
[thisc, clabel] = gui.i_selectnclass(sce,allowunique,[],[],FigureHandle);
if isempty(thisc), return; end
thisc = string(thisc);

% Shared with gui.callback_CompareCellScoreBtwCls so both group choosers
% look and behave the same - see gui.i_selectgroupsubset.
picked = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
if isempty(picked), return; end
if ~all(picked), thisc = thisc(picked); end

answer = gui.myQuestdlg(FigureHandle, "Violinplot for gene expression or cell state variables?","", ...
{'Gene Expression', 'Cell State'},'Gene Expression');
switch answer
    case 'Gene Expression'
        [glist] = gui.i_selectngenes(sce, [], FigureHandle);
        if isempty(glist)
            gui.myHelpdlg(FigureHandle, 'No gene selected.');
            return;
        end

        [Xt] = gui.i_transformx(sce.X, true, 9, FigureHandle);
        if isempty(Xt), return; end

        % if isscalar(glist)
        % answer='No';
        % else
        %    answer = gui.myQuestdlg(FigureHandle, 'Plot all in the same figure?','');
        % end

        % if strcmp(answer, 'Yes')
        %     [c, cL, noanswer] = gui.i_reordergroups(thisc, [], FigureHandle);
        %     if noanswer, return; end
        %     fw = gui.myWaitbar(FigureHandle);
        %
        %     gui.i_violinmatrix(full(Xt(:,picked)), sce.g, c, cL, glist, ...
        %             false, '', FigureHandle);
        %
        %     gui.myWaitbar(FigureHandle, fw);
        %
        %     return;
        % end

        n = length(glist);
        thisyv = cell(n,1);
        for k=1:n
            thisyv{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
            if ~all(picked)
                thisyv{k} = thisyv{k}(picked);
            end
        end
        ylabelv = glist;
        ylab = 'Expression';

    case 'Cell State'
        [thisyv, ylabelv] = gui.i_selectnstates(sce, true, [1], FigureHandle);

        a = false(length(thisyv), 1);
        for k = 1:length(thisyv)
            thisyv{k} = thisyv{k}(picked);
            a(k) = isnumeric(thisyv{k});
        end
        if any(a)
            if ~all(a)
                thisyv = thisyv(a);
                ylabelv = ylabelv(a);
                gui.myHelpdlg(FigureHandle, 'Only continuous variables of cell state will be shown.');
            end
        else
            gui.myHelpdlg(FigureHandle, 'No valid cell state variables. Violinplot cannot be shown.');
            return;
        end
        ylab = 'Value';
    otherwise
        return;
end
% assignin("base","ylabelv",ylabelv);
% assignin("base","thisyv",thisyv);
% assignin("base","thisc",thisc);
% The tab names the gene or state; these say what the axes are. Without the
% ylabel the plot defaults to 'Score', which is right for a cell score and
% wrong for expression.
labelinfo = struct('ylabel', ylab, 'xlabel', clabel);
gui.sc_uitabgrpfig_vioplot(thisyv, ylabelv, thisc, FigureHandle, [], labelinfo);
end
