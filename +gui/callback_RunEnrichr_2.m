function callback_RunEnrichr_2(src, ~, predefinedlist)

if nargin < 3, predefinedlist = []; end

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    gsorted = natsort(sce.g);

    rng("shuffle");
    n = length(gsorted);
    if isempty(predefinedlist)
        ingenelist = gui.i_inputgenelist(gsorted(randperm(n, min([200, length(gsorted)]))));
    else
        ingenelist = gui.i_inputgenelist(predefinedlist);
    end
    if isempty(ingenelist), return; end

    answer = questdlg('Add background list?','');
    switch answer
        case 'Yes'
            [idx] = gui.i_selmultidlg(sce.g, sce.g, FigureHandle);
            if isempty(idx), return; end
            if idx == 0, return; end
            backgroundlist = sce.g(idx);
        case 'No'
            backgroundlist = [];
        case 'Cancel'
            return;
    end

    [genesets] = in_selDataSources;
    if isempty(genesets), return; end

     definput = {'5', '0.1'};
    prompt = {'Min # of overlapping genes:', ...
              'P-value cutoff:'};
    dlgtitle = 'Enrichr Result Filter';
    dims = [1, 80];
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    if isempty(answer), return; end
    try
        minugenes = str2double(answer{1});
        pvaluecut = str2double(answer{2});
        assert((minugenes > 0) && (minugenes < 100));
        assert((pvaluecut >= 0.0) && (pvaluecut <= 1.0));
    catch
        errordlg('Invalid input.');
        return;
    end
    

    fw = gui.gui_waitbar;
    Tlist = run.mt_Enrichr(ingenelist, backgroundlist, genesets,...
                           minugenes, pvaluecut);

    T=table;
    for k = 1:height(Tlist)
        if ~isempty(Tlist{k}) && istable(Tlist{k})
            T = [T; Tlist{k}];
        end
    end

    gui.gui_waitbar(fw);
    
    [~, ~] = gui.i_exporttable(T, true, 'Tenrichrres', 'Enrichr_Results');
    gui.i_viewtable(T, FigureHandle);


    function [genesets] = in_selDataSources
        genesets = [];
        % dsv = ["GO_Biological_Process_2023", "GO_Molecular_Function_2023", ...
        %        "KEGG_2019_Mouse", "KEGG_2021_Human", "Reactome_2022", ... 
        %        "WikiPathways_2019_Mouse", "WikiPathways_2023_Human", ...
        %        "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019", ...
        %        "miRTarBase_2017"];
        dsv = pkg.i_get_enrichr_libraries;
        [idx1] = gui.i_selmultidlg(dsv, ["GO_Biological_Process_2023", ...
            "GO_Molecular_Function_2023"], FigureHandle);
        if isempty(idx1), return; end
        if idx1 == 0, return; end
        genesets = dsv(idx1);        
    end

end
