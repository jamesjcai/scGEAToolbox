function callback_FindAllMarkers(src, ~)

[~, sce, isui] = gui.gui_getfigsce(src);

    answer = questdlg('Select Method', ...
        '', 'Marker Gene Heatmap', 'Find All Markers', ...
        'Marker Gene Heatmap');
    switch answer
        case 'Find All Markers'
            in_findAllMarkers(sce);
        case 'Marker Gene Heatmap'
            in_MarkerGeneHeatmap(src);
    end
end


function in_findAllMarkers(sce)
    [thisc, ~] = gui.i_select1class(sce);
    if isempty(thisc), return; end
    if isscalar(unique(thisc))
        warndlg("All cells are in the same group.",'');
        return;
    end    
    [T] = pkg.e_findallmarkers(sce.X, sce.g, thisc, [], [], [], true);
    if ~isempty(T)
        needwait = true;
        [answer, filename] = gui.i_exporttable(T, needwait, 'Tallmarkers', ...
            'AllMarkersTable');
                % "Tcellattrib","CellAttribTable"
                % "Tviolindata","ViolinPlotTable"
                % "Tcrosstabul","CrosstabulTable"
                % "Tcellsignmt","CellSignatTable"
                % "Tdpgenesres","DPGenesResTable"
                % "Tallmarkers","AllMarkersTable"
        if ~isempty(answer)
            disp(filename);
            helpdlg(sprintf('All Markers Table saved.'), '');
        end
    else
        helpdlg('No results.', '');
    end
end

function in_MarkerGeneHeatmap(src, ~, sce)
    mfolder = fileparts(mfilename('fullpath'));
    if nargin < 3
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end
    % unique(sce.c_cluster_id)

    
    answer = questdlg("Only consider known (PangloaDB) marker genes?","");
    if strcmp(answer, 'Yes')
        markergenes = pkg.i_get_panglaodbmarkers;
        idx = ismember(upper(sce.g), upper(markergenes));
        if sum(idx) > 2000
            fprintf('Size of input matrix: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);
            sce.X = sce.X(idx, :);
            sce.g = sce.g(idx);            
            sce = sce.qcfilter;
            fprintf('Size of filtered matrix: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);
        else
            uiwait(helpdlg("Not enough marker genes (n < 2000). Use all genes to search for markers.",""));
        end
    elseif strcmp(answer, 'No')
        disp('Consider all genes for marker gene search.');
    else
        return;
    end
    
    [thisc, ~] = gui.i_select1class(sce);
    if isempty(thisc), return; end
    if isscalar(unique(thisc))
        warndlg("All cells are in the same group.",'');
        return;
    end
    % [c, cL, noanswer] = gui.i_reordergroups(thisc);
    % if noanswer, return; end
    
    [c] = grp2idx(thisc);
    answer = questdlg('Generate marker gene heatmap', ...
        'Select Method', 'Method 1 (DE 🐇)', 'Method 2 (scGeneFit 🐢)', ...
        'Method 3 (LASSO 🐢🐢)', 'Method 1 (DE 🐇)');
    switch answer
        case 'Method 1 (DE 🐇)'
            methodid = 1;
        case 'Method 2 (scGeneFit 🐢)'
            methodid = 2;
        case 'Method 3 (LASSO 🐢🐢)'
            methodid = 3;
        otherwise
            return;
    end
    
    fw = gui.gui_waitbar;

    % speciestag = gui.i_selectspecies(2, false, FigureHandle);
    % if isempty(speciestag)
    %     requirerefresh = false;
    %     return;
    % end
    speciestag = "human";
    load(fullfile(mfolder, ...
        '../resources', 'Biomart', sprintf('Biomart_%s_genes.mat',speciestag)), 'T');
    ApprovedSymbol = string(T.GeneName);
    sce = sce.rmmtgenes;
    sce = sce.rmhemoglobingenes;
    sce = sce.rmribosomalgenes;
    sce = sce.rmlncrnagenes;
    [idx] = ~ismember(upper(sce.g), upper(ApprovedSymbol));
    if any(idx)
        sce.g(idx) = [];
        sce.X(idx, :) = [];
    end
    
    fprintf('Size of used matrix*: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);

    try
        [markerlist] = sc_pickmarkers(sce.X, sce.g, c, 10, methodid);
    catch ME
        gui.gui_waitbar(fw, true);
        errordlg(ME.message);
        return;
    end
    
    glist = string(markerlist{1}(:));
    for k = 2:length(markerlist)
        glist = [glist; string(markerlist{k}(:))];
    end
    gui.gui_waitbar(fw);
    gui.i_heatmap(sce, glist, thisc, FigureHandle);    
end

