function xcallback_ShowPseudoTimeGenes(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[thisc, ~, ~, newpickclable] = gui.i_select1state(sce, true, false);
if isempty(thisc)
    helpdlg('No pseudotime values are available for cells.', '')
        return;
    end
    if ~isempty(newpickclable), clable = newpickclable; end

    fw = gui.gui_waitbar;
    try
        Xi = sc_impute(sce.X, 'type', 'MAGIC');
        [r, pval] = corr(thisc, Xi');
        [~, ~, ~, fdr] = pkg.fdr_bh(pval);
        T = table(sce.g, r', pval', fdr');
        T.Properties.VariableNames = {'genes', 'r', 'pval', 'fdr'};
        [T, idx] = sortrows(T, 'r', 'descend');
    catch ME
        gui.gui_waitbar(fw, true);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw, true);
    if ~(ismcc || isdeployed)
        labels = {'Save correlation table to variable named:'};
        vars = {'T'};
        values = {T};
        waitfor(export2wsdlg(labels, vars, values));
    else
        gui.i_exporttable(T, false, 'CORRTable');
    end

    answer = questdlg('Select genes to show correlation?', '');

    switch answer
        case 'Yes'
            sortedg = sce.g(idx);
            sortedX = Xi(idx, :);
            sortedgtxt = sortedg;
            for k = 1:length(sortedg)
                sortedgtxt(k) = sprintf('%s %.2f', T.genes(k), T.r(k));
            end

            [idx] = gui.i_selmultidlg(sortedgtxt);
            if isempty(idx), return; end
            if idx == 0, return; end
            glist = sortedg(idx);
            [c] = grp2idx(thisc);
            for k = 1:length(glist)
                f = figure('visible', 'off');
                scatter(c, sortedX(idx(k), :)');
                title(glist(k))
                xlabel('Pseudotime (stretched)')
                ylabel('Expression Level (imputated by MAGIC)')
                P = get(f, 'Position');
                set(f, 'Position', [P(1) - 20 * k, P(2) - 20 * k, P(3), P(4)]);
                set(f, 'visible', 'on');
                drawnow;
            end
    end
end
