function i_genescatter(T, parentfig)

    if nargin < 2, parentfig = []; end
    hx = gui.myFigure(parentfig);

    g = T.gene;

    hFig = hx.FigHandle;
    hAx = hx.AxHandle;
    if isempty(hAx), hAx = gca; end

    hx.addCustomButton('off', @in_callback_HighlightGenes, 'plotpicker-qqplot.gif', 'Highlight top HVGs');
    hx.addCustomButton('off', @in_HighlightSelectedGenes, 'curve-array.jpg', 'Highlight selected genes');
    hx.addCustomButton('off', @ExportGeneNames, 'bookmark-book.jpg', 'Export Selected HVG gene names...');
    hx.addCustomButton('off', @ExportTable, 'floppy-disk-arrow-in.jpg', 'Export HVG Table...');
    hx.addCustomButton('off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');
    hx.addCustomButton('off', @ChangeAlphaValue, 'Brightness-3--Streamline-Core.jpg', 'Change MarkerFaceAlpha value');
    
    h = scatter(hAx, T.de_coef, T.dv_coef, 'filled', 'MarkerFaceAlpha', .2);
    xlabel(hAx, 'DE Coefficient');
    ylabel(hAx, 'DV Coefficient');

    if ~isempty(g)
        dt = datacursormode(hFig);
        dt.UpdateFcn = {@in_myupdatefcn3, g};
    end
    hx.show;

    function ChangeAlphaValue(~, ~)
        if h.MarkerFaceAlpha <= 0.05
            h.MarkerFaceAlpha = 1;
        else
            h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
        end
    end

    function in_HighlightSelectedGenes(~,~)
        [glist] = gui.i_selectngenes(SingleCellExperiment([], T.gene),...
            [], hFig);
        if ~isempty(glist)            
            [y,idx]=ismember(glist, T.gene);
            idx=idx(y);            
            % idv = zeros(1, length(hvgidx));
            % idv(idx)=1;
            % h.BrushData = idv;
            for k=1:length(idx)
                dt = datatip(h,'DataIndex',idx(k));
            end
        end
    end

    function in_callback_HighlightGenes(~, ~)

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
            gui.myWarndlg(hFig, "No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n', sum(ptsSelected));

        gselected=g(ptsSelected);
        [yes,idx]=ismember(gselected, g);
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
            gui.myWarndlg(hFig,"No gene is selected.");
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
            txt = g(idx);
        else
            txt = num2str(event_obj.Position(2));
        end
    end
    
end