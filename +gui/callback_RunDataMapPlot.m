function callback_RunDataMapPlot(src, ~)

    if ~isa(src,'SingleCellExperiment')
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end
    if size(sce.s, 2) ~= 2
        gui.myWarndlg(FigureHandle, 'This function only works for 2D embedding.');
        return;
    end
    
    [thisc, ~] = gui.i_select1class(sce, true,[],[],FigureHandle);
    
    extprogname = 'py_datamapplot';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    
    try
        run.py_datamapplot(sce, thisc, wkdir);
    catch ME
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    end
    
    answer=gui.myQuestdlg(FigureHandle, ...
        sprintf('Script.py and input files saved. Open the folder %s?', wkdir));
    if strcmp(answer,'Yes'), winopen(outdir); end
end