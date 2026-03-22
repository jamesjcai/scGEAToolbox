function [requirerefresh] = callback_RunSeuratWorkflow(src, ~)

requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

extprogname = 'R_Seurat';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
[ok] = gui.i_confirmscript('Run Seurat Workflow (Seurat/R)?', ...
    'R_Seurat', 'r', FigureHandle);
if ~ok, return; end

[ndim] = gui.i_choose2d3d(FigureHandle);
if isempty(ndim), return; end
fw = gui.myWaitbar(FigureHandle);
try
    [sce] = run.r_seurat(sce, ndim, wkdir, true);
catch ME
    gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;
end
