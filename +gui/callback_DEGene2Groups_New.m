function callback_DEGene2Groups_New(src, ~)
    
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    % if ~gui.gui_showrefinfo('DE Analysis', FigureHandle), return; end

    extprogname = 'scgeatool_DEAnalysis';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    outdir = wkdir;

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


    outfile = sprintf("%s_vs_%s_DE_results.xlsx", ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)));

    filesaved = fullfile(outdir, outfile);
    try

        %[~, filesaved] = gui.i_exporttable(T, true, ...
        %    'Tdegenelist', outfile, 'Excel file', "All_raw", FigureHandle);        
        writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_genes');        
    catch
        return;
    end
    if ~isfile(filesaved), return; end
    % gui.myHelpdlg(FigureHandle, sprintf('Result has been saved in %s', filesaved));

    if ~strcmp('Yes', gui.myQuestdlg(FigureHandle, ...
            sprintf('Result has been saved in %s. Additional Analysis?', filesaved)))        
        if strcmp('Yes', gui.myQuestdlg(FigureHandle,'Open Output Folder?'))
            winopen(fileparts(filesaved));
        end
        return;
    end

    items = {'Set Filter Parameters', 'Enrichr Analysis', ...
        'LLM Summarize', 'Generate Volcano Plot', 'Open Output Folder'};
    selected = gui.myChecklistdlg(FigureHandle, items, ...
        'Title', 'Select Items','DefaultSelection', [2 4 5]);
    if isempty(selected)
        %writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_genes');
        %gui.myHelpdlg(FigureHandle, sprintf('Result has been saved in %s', filesaved));
        return;
    end


    % Process the selected analyses
    if any(contains(selected, 'Set Filter Parameters'))
        [paramset] = gui.i_degparamset(false, FigureHandle);
    else
        % paramset = [];
        preftagname ='degtestparamset';
        paramset = getpref('scgeatoolbox', preftagname, {0.05, 1.0, 0.01, 'Adjusted P-value'});
    end

    fw = gui.myWaitbar(FigureHandle);

    [Tup, Tdn, paramset] = pkg.e_processdetable(T, paramset, FigureHandle);
    [T, Tnt] = pkg.in_DETableProcess(T, cL1, cL2, sum(i1), sum(i2));

    try
        writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_processed');
        writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
        writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');
        writetable(Tnt, filesaved, "FileType", "spreadsheet", 'Sheet', 'Note');
    catch ME
        warning(ME.message);
    end

    
    % Perform additional analyses based on user selection
    if any(contains(selected, 'Enrichr Analysis'))
        gui.myWaitbar(FigureHandle, fw, false, '', 'Enrichr Analysis...', 2/3);
        try
            gui.e_enrichrxlsx(Tup,Tdn,T,filesaved);
        catch ME
            warning(ME.message);
        end
    end

    if any(contains(selected, 'LLM Summarize'))
        gui.myWaitbar(FigureHandle, fw, false, '', 'LLM Summarize...', 2/3);
        try
            [TbpUpEnrichr, TmfUpEnrichr, ...
                    TbpDnEnrichr, TmfDnEnrichr] = pkg.in_XLSX2DETable(filesaved);
            [~, wordfilename] = fileparts(filesaved);
        
            llm.e_DETableSummary(TbpUpEnrichr, ...
                TmfUpEnrichr, TbpDnEnrichr, ...
                TmfDnEnrichr, wordfilename, outdir);
        catch ME
            warning(ME.message);
        end
    end

    f = [];
    if any(contains(selected, 'Generate Volcano Plot'))
        gui.myWaitbar(FigureHandle, fw, false, '', 'Generate Volcano Plot...', 2/3);
        try
            f = e_volcano(T, Tup, Tdn, FigureHandle);
        catch ME
            warning(ME.message);
        end
    end   
    gui.myWaitbar(FigureHandle, fw);

    pause(2);

    if any(contains(selected, 'Open Output Folder'))        
        if isempty(f), f = FigureHandle; end
        if strcmp('Yes', gui.myQuestdlg(f,'Open Output Folder?'))
            winopen(fileparts(filesaved));
        end
    end

    function in_callback_savetable(srcx, ~)
        hFig = srcx.Parent.Parent;
        gui.i_exporttable(T, true, ...
            'Tdegenelist', outfile, [], "All_genes", hFig);
    end
   
    function in_callback_runenrichr(srcx, ~)
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
                outbackgroundlist, "Up", outdir);
        end
        
        [outgenelist, outbackgroundlist, enrichrtype] = ...
            gui.gui_prepenrichr(Tdn.gene, sce.g,... 
           'Run enrichment analysis with down-regulated DE genes?', ...
           hFig);
    
        if ~isempty(outbackgroundlist)
            gui.callback_RunEnrichr(src, [], outgenelist, enrichrtype, ...
                outbackgroundlist, "Down", outdir);
        end        
    end

% ------- begining of volcano_plot

    function hFig = e_volcano(T, Tup, Tdn, parentfig)
        T=T(~ismember(T.gene, [Tup.gene; Tdn.gene]),:);    
        hx = gui.myFigure(parentfig);
        hFig = hx.FigHandle;
        ax = hx.AxHandle;
        if isempty(ax), ax = gca; end
    
        hx.addCustomButton('off', @in_callback_runenrichr, 'www.jpg', 'Run enrichment analysis');
        hx.addCustomButton('off', @in_callback_savetable, 'floppy-disk-arrow-in.jpg', 'Export DE Gene Table...');
        
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

    % function i_replaceinf(T)
    %     % Iterate over each variable in the table
    %     for varIdx = 1:width(T)
    %         % Check if the variable is numeric
    %         if isnumeric(T{:, varIdx})
    %             % Replace positive Inf with 1e99
    %             T{T{:, varIdx} == Inf, varIdx} = 1e99;
    %             % Replace negative Inf with -1e99
    %             T{T{:, varIdx} == -Inf, varIdx} = -1e99;
    %         end
    %     end
    % 
    % end

end