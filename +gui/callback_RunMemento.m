function callback_RunMemento(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    % if ~gui.gui_showrefinfo('Memento [PMID:39454576]', FigureHandle), return; end

    %[wkdir] = gui.i_getwrkdir;
    %if isempty(wkdir), return; end
    extprogname = 'py_memento';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    
    [i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
    if isscalar(i1) || isscalar(i2), return; end
    
    % --------
    a=sprintf('%s vs. %s',cL1{1}, cL2{1});
    b=sprintf('%s vs. %s',cL2{1}, cL1{1});
    answer = gui.myQuestdlg(FigureHandle, 'Which vs. which?','',{a,b},a);
    switch answer
        case a
        case b
            i3=i1; i1=i2; i2=i3;
            cL3=cL1; cL1=cL2; cL2=cL3;
        otherwise
            return;
    end
    % ----------
    X1 = sce.X(:, i1);
    X2 = sce.X(:, i2);
    c = [zeros(size(X1,2),1); ones(size(X2,2),1)];
    scex = SingleCellExperiment([X1 X2], sce.g, [], c);
    scex.c_batch_id = c;
    [succeeded] = run.py_writeh5ad(scex, 'input.h5ad', wkdir);
    if ~succeeded, return; end


        T = run.py_memento(wkdir);
        outfile = sprintf('%s_vs_%s_Memento_results', ...
            matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
        [filetype, filesaved] = gui.i_exporttable(T, true, ...
            'Tmementores', outfile, [], "All_genes", FigureHandle);    
    
        if ~strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Generate scatter plot?')), return; end

    hx=gui.myFigure(FigureHandle);
    hFig = hx.FigHandle;
    hAx = hx.ax;

    hx.addCustomButton('off', @HighlightGenes, 'plotpicker-qqplot.gif', 'Highlight top HVGs');
    hx.addCustomButton('off', @in_HighlightSelectedGenes, 'curve-array.jpg', 'Highlight selected genes');
    hx.addCustomButton('off', @ExportGeneNames, 'bookmark-book.jpg', 'Export Selected HVG gene names...');
    hx.addCustomButton('off', @ExportTable, 'floppy-disk-arrow-in.jpg', 'Export HVG Table...');
    hx.addCustomButton('off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');
    hx.addCustomButton('off', @ChangeAlphaValue, 'Brightness-3--Streamline-Core.jpg', 'Change MarkerFaceAlpha value');
    
    h = scatter(hAx, T.de_coef, T.dv_coef, 'filled', 'MarkerFaceAlpha', .1);
    xlabel(hAx, 'DE Coefficient');
    ylabel(hAx, 'DV Coefficient');

    %if ~isempty(g)
        dt = datacursormode(hFig);
        dt.UpdateFcn = {@in_myupdatefcn3, string(T.gene)};
    %end
    hx.show;

    function ChangeAlphaValue(~, ~)
        if h.MarkerFaceAlpha <= 0.05
            h.MarkerFaceAlpha = 1;
        else
            h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
        end
    end

    function in_HighlightSelectedGenes(~,~)
        %Myc, Oct3/4, Sox2, Klf4
        [glist] = gui.i_selectngenes(SingleCellExperiment(X,g),...
            intersect(upper(g),["MYC", "POU5F1", "SOX2", "KLF4"]));
        if ~isempty(glist)            
            [y,idx]=ismember(glist,g);
            idx=idx(y);            
            % idv = zeros(1, length(hvgidx));
            % idv(idx)=1;
            % h.BrushData = idv;
            for k=1:length(idx)
                dt = datatip(h,'DataIndex',idx(k));
            end
        end
    end

    function HighlightGenes(~, ~)

        Tx = T(T.dv_coef > 0 & T.de_coef > 0,:);
        [~, hvgidx] = sort(Tx.dv_pval);
        
        idx = zeros(1, length(hvgidx));
        h.BrushData = idx;
        k = gui.i_inputnumk(200, 1, 2000);
        if isempty(k), return; end
        idx(hvgidx(1:k)) = 1;
        h.BrushData = idx;
    end

    function ExportTable(~, ~)
        gui.i_exporttable(T, true, 'Tmementores', 'MementoRsTable');
        % Tdegenelist 
        % 'Tviolindata','ViolinPlotTable'
        % 'Thvgreslist', 'HVGResultTable' 
    end

    function ExportGeneNames(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warning("No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n', sum(ptsSelected));


        gselected=g(ptsSelected);
        [yes,idx]=ismember(gselected,T.genes);
        Tx=T(idx,:);
        Tx=sortrows(Tx,4,'descend');
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.genes;


        labels = {'Save selected gene names to variable:',...
            'Save HVG table:'};
        vars = {'g','T'};
        values = {tgenes,T};
        export2wsdlg(labels, vars, values, ...
            'Save Data to Workspace');
    end

    function EnrichrHVGs(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warning("No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n', sum(ptsSelected));

        gselected=g(ptsSelected);
        [yes,idx]=ismember(gselected,T.genes);
        Tx=T(idx,:);
        Tx=sortrows(Tx,4,'descend');
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.genes;
        
        gui.i_enrichtest(tgenes, g, numel(tgenes));
    end
        
       
    
   
    function txt = in_myupdatefcn3(src, event_obj, g)  
    
        if isequal(get(src, 'Parent'), hAx)
            idx = event_obj.DataIndex;
            txt = {g(idx)};
        else
            txt = num2str(event_obj.Position(2));
        end
    end
    
    





end