function callback_DPGene2Groups(src, ~)

    [FigureHandle, sce_ori] = gui.gui_getfigsce(src);
    sce = copy(sce_ori);
    
    if ~gui.gui_showrefinfo('DP Analysis', FigureHandle), return; end
    
    extprogname = 'scgeatool_DPAnalysis';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end
    
    
    [i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
    if isscalar(i1) || isscalar(i2), return; end


    c=zeros(size(i1));
    c(i1)=1; c(i2)=2;
    cL=[cL1; cL2];
    
    if ~all(c>0)
        sce = sce.selectcells(c>0); % OK
        c = c(c>0);
        i1 = c==1;
        i2 = c==2;
    end
    
    [indx1,speciesid] = gui.i_selgenecollection(FigureHandle);
    if isempty(indx1), return; end
    [setmatrx, setnames, setgenes] = pkg.e_getgenesets(indx1, speciesid, FigureHandle); % (indx1);
    if isempty(setmatrx) || isempty(setnames) || isempty(setgenes)
        return;
    end
    
    % assignin("base", "setmatrx", setmatrx);
    % assignin("base", "setnames", setnames);
    % assignin("base", "setgenes", setgenes);
    
    fw = gui.myWaitbar(FigureHandle, [], false, 'Computing differential programs...');
    cleanupObj = onCleanup(@() i_closewaitbar(fw)); 
    
    ranknorm   = true;
    bgsubtract = true;
    try
        sceX = log1p(sc_norm(sce.X));
        % sceX = sc_transform(sce.X, 'type', 'PearsonResiduals');
        [T, Zx_dpg, Zy_dpg] = sc_dpg(sceX(:,i1), sceX(:,i2), sce.g, setmatrx, setnames, setgenes, ranknorm, bgsubtract);
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    
    [~, ix, iy] = intersect(upper(setgenes), upper(sce.g));
    setgenes = setgenes(ix);
    setmatrx = setmatrx(:,ix);
    sce.X = sce.X(iy,:);
    sce.g = sce.g(iy);
    
    gui.myWaitbar(FigureHandle, fw, false, '', 'Saving DP results...', 0.70);
    
    if height(T)==0
        i_closewaitbar(fw);
        gui.myHelpdlg(FigureHandle, 'No significant results.');
        return;
    end
    
    outfile = sprintf('%s_vs_%s_DP_results', ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)));
    
    gui.myWaitbar(FigureHandle, fw, false, '', 'Preparing DP results...', 0.85);
    i_closewaitbar(fw);
    
    plotAction = struct('Text', 'Plot Gene Sets', ...
        'Tooltip', 'Select rows in the table, then click to generate dot/violin plots', ...
        'Callback', @in_callback_plotgenesets);
    gui.TableViewerApp(T, FigureHandle, outfile, plotAction);


    function i_closewaitbar(fw)
        if nargin < 1 || isempty(fw)
            return;
        end
        
        try
            if pkg.i_isvalid(fw)
                close(fw);
            end
        catch
            % waitbar may already be closed; safe to ignore
        end
    end
    
    function in_callback_plotgenesets(tableObj, figtab)
        Tview = i_tableFromUITable(tableObj);
        if isempty(Tview) || ~ismember('setnames', Tview.Properties.VariableNames)
            gui.myWarndlg(figtab, 'No DP gene sets are available for plotting.');
            return;
        end
        
        if isempty(tableObj.Selection)
            gui.myWarndlg(figtab, 'Select one or more rows in the table first, then click Plot Gene Sets.');
            return;
        end
        selectedRows = unique(tableObj.Selection(:, 1));

        allSetNames = string(Tview.setnames);
        selectedSets = allSetNames(selectedRows);
        selectedSets = unique(selectedSets(strlength(selectedSets) > 0), 'stable');
        if isempty(selectedSets)
            gui.myWarndlg(figtab, 'None of the selected rows contain a gene set name.');
            return;
        end
        
        outdir = tempdir;
        Xt = sceX(iy, :);
        images = {};
        
        fw = gui.myWaitbar(figtab, [], false, 'Generating plots...');
        cleanupPlot = onCleanup(@() i_closewaitbar(fw)); %#ok<NASGU>
        allSaved = true;
        anySaved = false;
        
        for kk = 1:numel(selectedSets)
            setName = selectedSets(kk);
            idx = string(setnames) == setName;
            posg = string.empty(0, 1);
            posg = string(setgenes(setmatrx(idx, :)));
            posg = posg(strlength(posg) > 0);
            if isempty(posg)
                continue;
            end
        
            outfile1 = sprintf('dotplot_%s.png', matlab.lang.makeValidName(setName));
            outfile2 = sprintf('violplt_%s.png', matlab.lang.makeValidName(setName));
            filesaved1 = fullfile(outdir, outfile1);
            filesaved2 = fullfile(outdir, outfile2);
        
            frac = 0.75 + 0.24 * ((kk-1) ./ max(numel(selectedSets), 1));
            gui.myWaitbar(figtab, fw, false, '', ...
                sprintf('Generating plots (%d/%d)...', kk, numel(selectedSets)), frac);
        
            suc1 = false;
            try
                f1 = gui.i_dotplot(Xt, upper(sce.g), c, cL, upper(posg), true, setName);
                gui.i_exportfigure(f1, filesaved1);
                close(f1);
                images = [images, {filesaved1}];
                suc1 = true;
            catch ME
                warning(ME.message);
            end

            suc2 = false;
            try
                k_prog = find(string(setnames) == setName, 1);
                y = nan(numel(c), 1);
                y(i1) = Zx_dpg(k_prog, :)';
                y(i2) = Zy_dpg(k_prog, :)';
                f2 = gui.i_violinplot(y, cL(c), setName, true, [], posg, figtab);
                gui.i_exportfigure(f2, filesaved2);
                close(f2);
                images = [images, {filesaved2}];
                suc2 = true;
            catch ME
                warning(ME.message);
            end
        
            allSaved = allSaved && suc1 && suc2;
            anySaved = anySaved || suc1 || suc2;
        end
        
        gui.myWaitbar(figtab, fw, false, 'Finishing', 'Preparing PowerPoint export...', 0.99);
        if ~allSaved
            if anySaved
                gui.myHelpdlg(figtab, 'Some figure files could not be saved.');
            else
                gui.myHelpdlg(figtab, 'No figure files could be saved.');
            end
            winopen(outdir);
        end
        
        if ~isempty(images)
            gui.myWaitbar(figtab, fw, false, '', 'Exporting PowerPoint...', 0.995);
            gui.i_save2pptx(images, false, fw, figtab);
        end
    end
    
    function Tview = i_tableFromUITable(tableObj)
        columnNames = tableObj.ColumnName;
        data = tableObj.Data;
        Tview = cell2table(data, 'VariableNames', strrep(columnNames, ' ', '_'));
    end

end