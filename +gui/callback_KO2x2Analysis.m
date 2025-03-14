function callback_KO2x2Analysis(src, ~)
    if isa(src, 'SingleCellExperiment')
        sce = src;
        FigureHandle = [];
    else
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end
    % if ~gui.gui_showrefinfo('KO2x2 Analysis', FigureHandle), return; end

    extprogname = 'scgeatool_KO2x2Analysis';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    
    [thisc] = gui.i_select1class(sce, false, '','', FigureHandle);
    [c, cL] = grp2idx(thisc);
    cL = strrep(cL,'+','_pos');
    cL = strrep(cL,'-','_neg');
    if numel(unique(c)) ~= 4
        gui.myErrordlg(FigureHandle, 'Need exactly 4 categories of cells for this analysis.');
        return;
    end
    assert(length(cL)==4);

    fw = gui.myWaitbar(FigureHandle);
    count=0;
    for k = 1:4
        for l = k+1:4
            count=count+1;
            fprintf('Working on ...... %d\n', count);
            cL1=cL{k}; cL2=cL{l};
            gui.myWaitbar(FigureHandle, fw, false, '', ...
                sprintf('%s vs. %s', cL1, cL2), count/6);
            T = sc_deg(sce.X(:, c==k), sce.X(:, c==l), sce.g, 1, ...
                    false, FigureHandle);            
            in_processT(T, cL1, cL2);            
        end
    end
    gui.myWaitbar(FigureHandle, fw);


    function in_processT(T, cL1, cL2)
        try
            T = sortrows(T, 'p_val_adj', 'ascend');
            T = sortrows(T, 'pct_1', 'ascend');
            T = sortrows(T, 'pct_2', 'descend');
            T = sortrows(T, 'avg_log2FC', 'ascend');
        catch ME
            warning(ME.message);
        end
    
        try
            % avg_1 = mean(X,2);
            % avg_2 = mean(Y,2);
            % pct_1 = sum(X>0,2)./size(X,2);
            % pct_2 = sum(Y>0,2)./size(Y,2);
            if contains(T.Properties.VariableNames{5}, 'avg_1')
                T.Properties.VariableNames{5} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{5}, ...
                    matlab.lang.makeValidName(string(cL1)));
            end
    
            if contains(T.Properties.VariableNames{6}, 'avg_2')
                T.Properties.VariableNames{6} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{6}, ...
                    matlab.lang.makeValidName(string(cL2)));
            end
    
            if contains(T.Properties.VariableNames{7}, 'pct_1')
                T.Properties.VariableNames{7} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{7}, ...
                    matlab.lang.makeValidName(string(cL1)));
            end
    
            if contains(T.Properties.VariableNames{8}, 'pct_2')
                T.Properties.VariableNames{8} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{8}, ...
                    matlab.lang.makeValidName(string(cL2)));
            end
        catch ME
            warning(ME.message);
        end
    
        outfile = sprintf('%s_vs_%s_KO2x2_results', ...
            matlab.lang.makeValidName(string(cL1)), ...
            matlab.lang.makeValidName(string(cL2)));

        filesaved = fullfile(wkdir, outfile+".xlsx");
        
        writetable(T, filesaved, 'FileType', 'spreadsheet', ...
                    'WriteRowNames', true);
        
        paramset = {0.005, 0.1, 0.01, 'Adjusted P-value'};
        [Tup, Tdn] = pkg.e_processDETable(T, paramset, FigureHandle);
        if ~isempty(Tup)
            writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
        end
        if ~isempty(Tdn)
            writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
        end
        if ~isempty(Tup) || ~isempty(Tdn)
            gui.e_enrichrxlsx(Tup,Tdn,T,filesaved);
        end
        fprintf('\nFile saved: %s\n', filesaved);
    end

end
