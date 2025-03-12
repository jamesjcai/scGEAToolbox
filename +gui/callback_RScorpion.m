function callback_RScorpion(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('SCORPION [PMID:38438786]', FigureHandle), return; end

extprogname = 'R_SCORPION';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);

[ok] = gui.i_confirmscript('Using SCORPION to construct a GRN', ...
    extprogname, 'r');
if ~ok, return; end


fw = gui.myWaitbar(FigureHandle);
try
    T = run.r_SCORPION(sce.X,sce.g,wkdir,false);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
end

