function [needupdatesce] = callback_RunCopyKat(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

needupdatesce = false;
extprogname = 'R_copykat';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end


        [ok] = gui.i_confirmscript('Run R_copykat?', ...
            'R_copykat', 'r');
        if ~ok, return; end

        [ndim] = gui.i_choose2d3d;
        if isempty(ndim), return; end
        fw = gui.myWaitbar(FigureHandle);
        try
            run.r_copykat(sce, wkdir);
        catch
            gui.myWaitbar(FigureHandle, fw);
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
        guidata(FigureHandle, sce);
end