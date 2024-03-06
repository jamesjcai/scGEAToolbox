function callback_RScorpion(src, ~)

FigureHandle = src.Parent.Parent;
if ~gui.gui_uishowrefinfo('SCORPION [PMID:38438786]',FigureHandle)
    return;
end

extprogname = 'R_SCORPION';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);

[ok] = gui.i_confirmscript('Using SCORPION to construct a GRN', ...
    extprogname, 'r');
if ~ok, return; end


sce = guidata(FigureHandle);

fw = gui.gui_waitbar;
try
    T = run.r_SCORPION(sce.X,sce.g,wkdir,false);    
catch ME
    gui.gui_waitbar(fw, true);
    errordlg(ME.message);
    return;
end
gui.gui_waitbar(fw);
end

