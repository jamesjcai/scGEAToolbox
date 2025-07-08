function callback_DEGene2Groups(src, ~)
    
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    % if ~gui.gui_showrefinfo('DE Analysis', FigureHandle), return; end

    extprogname = 'SCGEATOOL_DEAnalysis';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    
    [i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
    if isscalar(i1) || isscalar(i2), return; end
    
    % --------
    a = sprintf('%s vs. %s',cL1{1}, cL2{1});
    b = sprintf('%s vs. %s',cL2{1}, cL1{1});
    answer = gui.myQuestdlg(FigureHandle, 'Which vs. which?','',{a,b}, a);
    switch answer
        case a
        case b
            i3=i1; i1=i2; i2=i3;
            cL3=cL1; cL1=cL2; cL2=cL3;
        otherwise
            return;
    end
    % ----------
    
    methodtag = "ranksum";   
    try
        switch methodtag
            case 'ranksum'
                T = sc_deg(sce.X(:, i1), sce.X(:, i2), sce.g, 1, ...
                    true, FigureHandle);
            case 'deseq2'
                [ok] = gui.i_confirmscript('DE analysis (DESeq2)', ...
                    'R_DESeq2', 'r');
                if ~ok, return; end
                fw = gui.myWaitbar(FigureHandle);
                T = run.r_DESeq2(sce.X(:, i1), sce.X(:, i2), sce.g);
                gui.myWaitbar(FigureHandle, fw);
            case 'mast'
                [ok] = gui.i_confirmscript('DE analysis (MAST)', 'R_MAST', 'r');
                if ~ok, return; end
                fw = gui.myWaitbar(FigureHandle);
                T = run.r_MAST(sce.X(:, i1), sce.X(:, i2), sce.g);
                gui.myWaitbar(FigureHandle, fw);
        end
    catch ME
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end   

    [T] = in_Tprocessing(T);

    outfile = sprintf("%s_vs_%s_DE_results.xlsx", ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)));
    
    
    didit = false;
    try
        tempfilesaved = fullfile(wkdir, outfile+".xlsx");
        writetable(i_replaceinf(T), tempfilesaved, "FileType", ...
            "spreadsheet", 'Sheet', 'All_genes');
        didit = true;
    catch ME
        warning(ME.message);
    end
    if ~didit
        try
            tempfilesaved = fullfile(wkdir, outfile+".csv");
            writetable(T, tempfilesaved);
        catch ME
            warning(ME.message);
        end
    end

    % gui.DEResultViewApp(T);
    [filetype, filesaved] = gui.i_exporttable(T, true, ...
        'Tdegenelist', outfile, [], "All_genes", FigureHandle);

    tf = 0;
    if ~(ismcc || isdeployed) && strcmp(filetype, 'Workspace')
        [Tup, Tdn, paramset] = pkg.e_processdetable(T, [], FigureHandle);
        if isempty(Tup) && isempty(Tdn), return; end
        labels = {'Save DE results (selected up-regulated) to variable named:', ...
            'Save DE results (selected down-regulated) to variable named:'};
        vars = {'Tup', 'Tdn'};
        values = {Tup, Tdn};
        [~, tf] = export2wsdlg(labels, vars, values);
        if tf ~= 1, return; end
    end

    if ~isempty(filesaved)
        if strcmp(filetype, 'Excel file')
            %answer = gui.myQuestdlg(FigureHandle, 'Save up- and down-regulated genes to seperate sheets?');
            %if strcmp(answer, 'Yes')
                [Tup, Tdn, paramset] = pkg.e_processdetable(T,[],FigureHandle);
                % strcmp(extractAfter(filesaved,strlength(filesaved)-4),'xlsx')
                if ~isempty(Tup)
                    writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
                end
                if ~isempty(Tdn)
                    writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
                end
                gui.myHelpdlg(FigureHandle, sprintf('Result has been saved in %s', filesaved));
                %writetable(Tup,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet',);
                %writetable(Tdn,fullfile(tempdir,sprintf('%s_up.xlsx',outfile)),'FileType','spreadsheet');
                tf = 1;
            %end
        elseif strcmp(filetype, 'Text file')
            [filePath, name, ext] = fileparts(filesaved);
            % strcmp(extractAfter(filesaved,strlength(filesaved)-3),'txt')
            answer = gui.myQuestdlg(FigureHandle, 'Save up- and down-regulated genes to seperate files?');
            if strcmp(answer, 'Yes')
                [Tup, Tdn, paramset] = pkg.e_processdetable(T,[],FigureHandle);
                if ~isempty(Tup)
                    namex = sprintf('%s_Upregulated', name);
                    f = fullfile(filePath, [namex, ext]);
                    if exist(f,"file"), f = 'Upregulated'; end
                    
                    [~, ~] = gui.i_exporttable(Tup, true, 'Tup', ...
                        f, 'Text file','', FigureHandle);
                end
                if ~isempty(Tdn)
                    namex = sprintf('%s_Downregulated', name);
                    f = fullfile(filePath, [namex, ext]);
                    if exist(f,"file"), f = 'Downregulated'; end
                                        
                    [~, ~] = gui.i_exporttable(Tdn, true, 'Tdn', ...
                        f, 'Text file','', FigureHandle);
                end
                tf = 1;
            end
        end
    end

    if tf ~= 1, return; end


    if contains(tempfilesaved, 'xlsx')
        if ~isempty(Tup)
            writetable(Tup, tempfilesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
        end
        if ~isempty(Tdn)
            writetable(Tdn, tempfilesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
        end    
    else
        if ~isempty(Tup)
            filename = fullfile(wkdir, outfile+"_Upregulated.csv");
            writetable(Tup, filename, 'WriteRowNames', true);            
        end
        if ~isempty(Tdn)
            filename = fullfile(wkdir, outfile+"_Downregulated.csv");
            writetable(Tdn, filename, 'WriteRowNames', true);
        end
    end

    if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Generate volcano plot?'))
        e_volcano(T, Tup, Tdn, FigureHandle);
        % exportgraphics(f, fullfile(wkdir, outfile+".png"), 'Resolution', 300);
    end

    function e_savetable(srcx, ~)
        hFig = srcx.Parent.Parent;
        [filetype, filesaved] = gui.i_exporttable(T, true, ...
            'Tdegenelist', outfile, [], "All_genes", hFig);
    end

    function [T] = in_Tprocessing(T)
        try
            T = sortrows(T, 'p_val_adj', 'ascend');
            T = sortrows(T, 'pct_1', 'ascend');
            T = sortrows(T, 'pct_2', 'descend');
            T = sortrows(T, 'avg_log2FC', 'ascend');
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
    end
    
    function e_runenrichr(srcx, ~)
        hFig = srcx.Parent.Parent;
        disp('To run enrichment analysis, type:');
        disp('run.web_Enrichr(Tup.gene(1:250))');
        disp('run.web_Enrichr(Tdn.gene(1:250))');
    
        [outgenelist, outbackgroundlist, enrichrtype] = ...
            gui.gui_prepenrichr(Tup.gene, sce.g,... 
           'Run enrichment analysis with up-regulated DE genes?', ...
           hFig);
    
        if ~isempty(outbackgroundlist)
            gui.callback_RunEnrichr(src, [], outgenelist, ...
                enrichrtype, ...
                outbackgroundlist, "Up", wkdir);
        end
        
        [outgenelist, outbackgroundlist, enrichrtype] = ...
            gui.gui_prepenrichr(Tdn.gene, sce.g,... 
           'Run enrichment analysis with down-regulated DE genes?', ...
           hFig);
    
        if ~isempty(outbackgroundlist)
            gui.callback_RunEnrichr(src, [], outgenelist, enrichrtype, ...
                outbackgroundlist, "Down", wkdir);
        end        
    end

% ------- begining of volcano_plot

    function hFig = e_volcano(T, Tup, Tdn, parentfig)
        T=T(~ismember(T.gene, [Tup.gene; Tdn.gene]),:);    
        hx = gui.myFigure(parentfig);
        hFig = hx.FigHandle;
        ax = hx.AxHandle;
        if isempty(ax), ax = gca; end
    
        hx.addCustomButton('off', @e_runenrichr, 'www.jpg', 'Run enrichment analysis');
        hx.addCustomButton('off', @e_savetable, 'floppy-disk-arrow-in.jpg', 'Export DE Gene Table...');
        
        hold(ax, "on");
        e_v(Tdn, ax);
        h = e_v(T, ax);
        e_v(Tup, ax);

        h.MarkerFaceColor=[.5 .5 .5];
        h.MarkerEdgeColor=[.5 .5 .5];
        title(ax, sprintf('%s vs. %s', ...
            string(cL1), ...
            string(cL2)));
        ylabel(ax, '-log_{10}(Adj. P-value)')
        xlabel(ax, 'log_{2}(FC)');
        lgd = legend(ax, {sprintf('Down-regulated (%d)', height(Tdn)), ...
            sprintf('Not Sig. (%d)', height(T)), ...
            sprintf('Up-regulated (%d)', height(Tup))},'Location', ...
            'bestoutside');
        try
            mindiffpct = paramset{1};
            minabsolfc = paramset{2};
            apvaluecut = paramset{3};
    
            Text_below_legend = sprintf('Dropout Diff. > %d%%\nLog2(FC) > %.2f\nAdj. P-Value < %g', ...
                100*mindiffpct, minabsolfc, apvaluecut);
            txt = text(ax, 0, 0, Text_below_legend, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', ...
                'FontSize', 10);
            % drawnow;
            updateTextbox;
            hFig.AutoResizeChildren = 'off';
            hFig.SizeChangedFcn = @updateTextbox;
        catch ME
            disp(ME.message);
        end
        hold(ax, "off");
        hx.show;
    
        function updateTextbox(~, ~)
            ax.Units = 'normalized';
            axPos = ax.Position;
            lgd.Units = 'normalized';
            lgdPos = lgd.Position;        
            legendCenterNorm = [lgdPos(1) + lgdPos(3)/2, ...
                                lgdPos(2)];
            axesNormX = (legendCenterNorm(1) - axPos(1)) / axPos(3);
            axesNormY = (legendCenterNorm(2) - axPos(2)) / axPos(4);
            xLimits = xlim(ax);
            yLimits = ylim(ax);        
            xData = xLimits(1) + axesNormX * (xLimits(2) - xLimits(1));            
            yData = yLimits(1) + axesNormY * (yLimits(2) - yLimits(1)) - 0.05 * range(yLimits);
            txt.Position = [xData, yData, 0];
        end

        function h = e_v(T, ax)
            genelist = T.gene; 
            pvals = T.p_val_adj;
            fc = T.avg_log2FC;
            h = ix_volcanoplot(fc, pvals, genelist, ax);
        end
        
        function h = ix_volcanoplot(fc, pvals, genelist, ax)
            %Vocano plot            
            pvals(pvals < 1e-100) = 1e-100;
            fc(fc<-999) = -10;
            fc(fc>999) = 10;
            x = fc;
            y = -log10(pvals);
            % [~, idx] = maxk(abs(y), 5);
            h = scatter(ax, x, y, 8, "filled");
            % hold(ax,"on");
            % scatter(x(idx),y(idx),'rx');
            %for k = 1:length(idx)
                %text(x(idx(k))+0.05, y(idx(k)), genelist(idx(k)));
            %end

            if isempty(genelist)
                disp('Empty genelist.');
            else
                h.DataTipTemplate.DataTipRows = dataTipTextRow('', genelist);
            end
            if any(abs(xlim(ax))>=10)
                xlim(ax,[-10 10]);
            end
        end
    end

    % ------- end of volcano_plot

    function [T] = i_replaceinf(T)
        % Iterate over each variable in the table
        for varIdx = 1:width(T)
            % Check if the variable is numeric
            if isnumeric(T{:, varIdx})
                % Replace positive Inf with 1e99
                T{T{:, varIdx} == Inf, varIdx} = 1e99;
                % Replace negative Inf with -1e99
                T{T{:, varIdx} == -Inf, varIdx} = -1e99;
            end
        end
        
    end

end