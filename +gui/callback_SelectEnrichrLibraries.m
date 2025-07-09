function callback_SelectEnrichrLibraries(src, ~)

    [FigureHandle] = gui.gui_getfigsce(src);
    try
        [genesets] = in_selDataSources;
    catch ME
        gui.myErrordlg(FigureHandle,'Runtime error');
    end
    if ~isempty(genesets)
        gui.myHelpdlg(FigureHandle, 'Enrichr gene set libraries selected.');
    end


    function [genesets] = in_selDataSources
        genesets = [];
        fw = gui.myWaitbar(FigureHandle);
        try
        
        dsv = pkg.i_get_enrichr_libraries;

        enrichrlibraries = getpref('scgeatoolbox', 'enrichrlibraries', ...
                                   ["GO_Biological_Process_2025", ...
                                    "GO_Molecular_Function_2025", ...
                                     "KEGG_2021_Human",...
                                     "Reactome_Pathways_2024"]);
        catch ME
            gui.myWaitbar(FigureHandle,fw,true);
            rethrow(ME);
        end
        gui.myWaitbar(FigureHandle,fw);
        [idx1] = gui.i_selmultidlg(dsv, enrichrlibraries, FigureHandle);
        if isempty(idx1), return; end
        if idx1 == 0, return; end
        genesets = dsv(idx1);
        setpref('scgeatoolbox', 'enrichrlibraries', genesets);    
    end

end
