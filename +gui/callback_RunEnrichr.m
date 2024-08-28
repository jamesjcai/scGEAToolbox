function callback_RunEnrichr(src, ~, predefinedlist)

if nargin < 3, predefinedlist = []; end

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    gsorted = natsort(sce.g);

    rng("shuffle");
    n = length(gsorted);
    if isempty(predefinedlist)
        tg = gui.i_inputgenelist(gsorted(randperm(n, min([200, length(gsorted)]))));
    else
        tg = gui.i_inputgenelist(predefinedlist);
    end
    if isempty(tg), return; end

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
    % if ~iscell(genesets), genesets = {genesets}; end
    fw = gui.gui_waitbar;

    Tlist = run.mt_Enrichr(tg, backgroundlist, genesets, 5);

    T=table;
    for k = 1:length(size(Tlist, 1))
        if ~isempty(Tlist{k, 1}) && istable(Tlist{k, 1})
            T = [T; Tlist{k, 1}];
        end
    end

    gui.gui_waitbar(fw);
    
    [filetype, filesaved] = gui.i_exporttable(T, true, 'Tenrichrres', 'Enrichr_Results');
    gui.i_viewtable(T, FigureHandle);


    function [genesets] = in_selDataSources
        genesets = [];
        dsv = ["GO_Biological_Process_2023", "GO_Molecular_Function_2023", ...
               "KEGG_2019_Mouse", "KEGG_2021_Human", "Reactome_2022", ... 
               "WikiPathways_2019_Mouse", "WikiPathways_2023_Human", ...
               "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019", ...
               "miRTarBase_2017"];
        [idx] = gui.i_selmultidlg(dsv, ["GO_Biological_Process_2023", ...
            "GO_Molecular_Function_2023"], FigureHandle);
        if isempty(idx), return; end
        if idx == 0, return; end
        genesets = dsv(idx);        
    end

end
