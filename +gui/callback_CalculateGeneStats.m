function callback_CalculateGeneStats(src, ~)
answer = questdlg('Calculate gene expression mean, CV, and dropout rate. Save output to a table.', '');
if ~strcmp(answer, 'Yes'), return; end

[~, sce] = gui.gui_getfigsce(src);

[Xt] = gui.i_transformx(sce.X);
if isempty(Xt), return; end


answer = questdlg('Grouping cells? Select Yes to pick a grouping variable. Select No to include all cells.', '');

if strcmp(answer, 'Yes')
    [thisc, ~] = gui.i_select1class(sce);
    if isempty(thisc), return; end
    [c, cL] = grp2idx(thisc);

    [newidx] = gui.i_selmultidlg(cL, natsort(cL));
    if isempty(newidx), return; end
    cx = c;
    c = zeros(size(c));
    for k = 1:length(newidx)
        c(cx == newidx(k)) = k;
    end
    cL = cL(newidx);
    
    fw = gui.gui_waitbar;
    T = sc_genestats(Xt(:, c == 1), sce.g);
    for j = 2:4
        T.Properties.VariableNames{j} = sprintf('%s_%s', T.Properties.VariableNames{j}, cL{1});
    end
    for k = 2:length(cL)
        [t] = sc_genestats(Xt(:, c == k), sce.g);
        for j = 2:4
            t.Properties.VariableNames{j} = sprintf('%s_%s', t.Properties.VariableNames{j}, cL{k});
        end
        T = [T, t(:, 2:4)];
    end
    gui.gui_waitbar(fw);
else
    fw = gui.gui_waitbar;
    T = sc_genestats(Xt, sce.g);
    gui.gui_waitbar(fw);
end

gui.i_exporttable(T,true,"Tgenestats","GeneStatsTable");
end