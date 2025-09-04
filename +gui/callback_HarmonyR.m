function [done] = callback_HarmonyR(src, ~)

done = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);
if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, 'No batch effect (all cells have the same SCE.C_BATCH_ID)');
    return;
end

extprogname = 'R_harmony';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

    %try
    id = sce.c_batch_id;
    if ~isnumeric(id)
        id = findgroups(sce.c_batch_id);
        id = id(:);
    end
    
    % [s] = run.py_harmonypy(sce.s, id, wkdir);
    [s] = run.r_harmony(sce.s, id, wkdir);

    if isempty(s) || isequal(sce.s, s)
        % gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, "Harmony Running Error");
        return;
    end
    sce.s = s;
    % catch ME
    %     %gui.myWaitbar(FigureHandle, fw,true);
    %     gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    %     %rethrow(ME);
    %     return;
    % end
    %gui.myWaitbar(FigureHandle, fw);
    gui.myGuidata(FigureHandle, sce, src);
    done = true;

end