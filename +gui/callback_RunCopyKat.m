function [needupdatesce] = callback_RunCopyKat(src, ~)

needupdatesce = false;
extprogname = 'R_copykat';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wkdir), return; end


        [ok] = gui.i_confirmscript('Run R_copykat?', ...
            'R_copykat', 'r');
        if ~ok, return; end

        [ndim] = gui.i_choose2d3d;
        if isempty(ndim), return; end
        fw = gui.gui_waitbar;
        try
            run.r_copykat(sce, wkdir);
        catch
            gui.gui_waitbar(fw);
            return;
        end
        gui.gui_waitbar(fw);
        guidata(FigureHandle, sce);
end