function [done] = callback_Harmonypy(src, ~)
%CALLBACK_HARMONYPY Batch integration with Harmony, Python backend.
%   Runs harmonypy, falling back to the native MATLAB implementation when
%   Python is not configured.
%
%   The backend question used to be asked here. It moved to
%   GUI.CALLBACK_HARMONYANY, the single menu entry that now reaches all
%   three backends, so reaching this file means Python was already chosen.
%
%   See also GUI.CALLBACK_HARMONYANY, GUI.CALLBACK_HARMONY, RUN.PY_HARMONYPY.

done = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);
if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, ...
        'No batch effect (all cells have the same SCE.C_BATCH_ID)');
    return;
end

if ~pkg.i_checkpython
    done = gui.callback_Harmony(src);
    return;
end

extprogname = 'py_harmonypy';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
if ~gui.i_setpyenv([], [], FigureHandle), return; end

id = sce.c_batch_id;
if ~isnumeric(id)
    id = findgroups(sce.c_batch_id);
end
id = id(:);

fw = gui.myWaitbar(FigureHandle);
try
    [s] = run.py_harmonypy(sce.s, id, wkdir);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(s) || isequal(sce.s, s)
    gui.myErrordlg(FigureHandle, "Harmonypy Running Error");
    return;
end

sce.s = s;
gui.myGuidata(FigureHandle, sce, src);
done = true;

end
