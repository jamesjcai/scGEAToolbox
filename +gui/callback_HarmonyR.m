function [done] = callback_HarmonyR(src, ~)
%CALLBACK_HARMONYR Batch integration with Harmony, R backend.
%   Offers the R package when R is configured, and the native MATLAB
%   implementation otherwise, so this menu entry works without R installed.
%
%   See also GUI.CALLBACK_HARMONY, RUN.R_HARMONY, RUN.ML_HARMONY.

done = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);
if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, ...
        'No batch effect (all cells have the same SCE.C_BATCH_ID)');
    return;
end

hasR = ispref('scgeatoolbox', 'rexecutablepath') && ...
    ~isempty(getpref('scgeatoolbox', 'rexecutablepath', []));
if ~hasR
    done = gui.callback_Harmony(src);
    return;
end

backend = gui.myQuestdlg(FigureHandle, 'Choose Harmony backend:', '', ...
    {'R (harmony package)', 'MATLAB (native)'}, 'R (harmony package)');
if isempty(backend), return; end
if strcmp(backend, 'MATLAB (native)')
    done = gui.callback_Harmony(src);
    return;
end

extprogname = 'R_harmony';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

id = sce.c_batch_id;
if ~isnumeric(id)
    id = findgroups(sce.c_batch_id);
end
id = id(:);

fw = gui.myWaitbar(FigureHandle);
try
    [s] = run.r_harmony(sce.s, id, wkdir);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(s) || isequal(sce.s, s)
    gui.myErrordlg(FigureHandle, "Harmony Running Error");
    return;
end

sce.s = s;
gui.myGuidata(FigureHandle, sce, src);
done = true;

end
