function [done] = callback_Harmonypy(src, ~)
%CALLBACK_HARMONYPY Batch integration with Harmony, Python backend.
%   Offers harmonypy when Python is configured, and the native MATLAB
%   implementation otherwise, so this menu entry works without Python
%   installed.
%
%   See also GUI.CALLBACK_HARMONY, RUN.PY_HARMONYPY, RUN.ML_HARMONY.

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

backend = gui.myQuestdlg(FigureHandle, 'Choose Harmony backend:', '', ...
    {'Python (harmonypy)', 'MATLAB (native)'}, 'Python (harmonypy)');
if isempty(backend), return; end
if strcmp(backend, 'MATLAB (native)')
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
