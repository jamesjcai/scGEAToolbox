function callback_DecontX(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~gui.gui_showrefinfo('DecontX [PMID:32138770]', FigureHandle)
    return;
end

extprogname = 'R_decontX';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
fw = gui.myWaitbar(FigureHandle);
try
    [Xdecon, ~] = run.r_decontX(sce, wkdir);
catch ME
    gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if ~(size(Xdecon, 1) == sce.NumGenes && ...
        size(Xdecon, 2) == sce.NumCells)
    gui.myErrordlg(FigureHandle, 'DecontX runtime error.','');
    return;
end

if strcmp(gui.myQuestdlg(FigureHandle, ...
        "Remove contamination? Click ''No'' " + ...
        "to keep data unchanged."), 'Yes')
    sce.X = round(Xdecon);
    gui.myGuidata(FigureHandle, sce, src);
    gui.myHelpdlg(FigureHandle, 'Contamination removed.');
end
end
