function callback_CalculateGeneStats(src, ~)
    % Callback function for Calculate Gene Stats
    % Prompts the user to calculate gene expression statistics and export the results.

    % Prompt user for calculation confirmation
    [FigureHandle, sce] = gui.gui_getfigsce(src);

    answer = gui.myQuestdlg(FigureHandle, 'Calculate gene expression mean, CV, and dropout rate. Save output to a table.', 'Confirmation');
    if ~strcmp(answer, 'Yes'), return; end

    % Retrieve single-cell experiment data
    Xt = gui.i_transformx(sce.X, [], [], FigureHandle);
    if isempty(Xt), return; end

    % Grouping cells based on user input
    groupingConfirmed = gui.myQuestdlg(FigureHandle, 'Grouping cells? Select Yes to pick a grouping variable. Select No to include all cells.', 'Grouping');
    if strcmp(groupingConfirmed, 'Yes')
        thisc = gui.i_select1class(sce);
        if isempty(thisc), return; end
        [c, cL] = grp2idx(thisc);

        newidx = gui.i_selmultidlg(cL, natsort(cL), FigureHandle);
        if isempty(newidx), return; end

        cx = c;
        c = zeros(size(c));
        for k = 1:length(newidx)
            c(cx == newidx(k)) = k;
        end
        cL = cL(newidx);
        
        % Initialize the table with default statistics
        T = sc_genestats(Xt(:, c == 1), sce.g);
        for j = 2:4
            T.Properties.VariableNames{j} = sprintf('%s_%s', T.Properties.VariableNames{j}, cL{1});
        end

        % Calculate and merge gene stats for each group
        for k = 2:length(cL)
            t = sc_genestats(Xt(:, c == k), sce.g);
            for j = 2:4
                t.Properties.VariableNames{j} = sprintf('%s_%s', t.Properties.VariableNames{j}, cL{k});
            end
            T = [T; t(:, 2:4)];
        end

    else
        % Calculate statistics for all cells if no grouping
        T = sc_genestats(Xt, sce.g);
    end
    
    % Export the results to a table
    gui.i_exporttable(T, true, 'GeneStatsTable',[],[],[], FigureHandle);
    
    % Update the GUI with waitbar
    gui.gui_waitbar;
end

