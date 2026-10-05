function [needupdatesce, T] = callback_RunPanhumanpy(src, ~)

needupdatesce = false;
T = [];   % the app always asks for T, and the early returns set none
[y, prepare_input_only] = gui.i_memorychecked([], []);
if ~y, return; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

% Preparing input files does not touch sce.c_cell_type_tx; only ask about
% overwriting labels when this run will actually assign new ones.
if ~prepare_input_only && ~gui.i_confirmoverwritecelltype(FigureHandle, sce)
    return;
end

% https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html
% SCimilarity trained model. Download SCimilarity models.
% Note, this is a large tarball - downloading and uncompressing can take a several minutes.

extprogname = 'py_panhumanpy';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

% panhumanpy matches upper-case symbols, so the run needs them, but SCE is
% the app's live handle: upper-casing it in place renamed every gene for
% the rest of the session (mouse Actb became ACTB), on the prepare-only
% path too. Restored on any return, as GUI.CALLBACK_RUNSCIMILARITY does.
originalGeneList = sce.g;
restoreGeneList = onCleanup(@() i_restoregenelist(sce, originalGeneList));
sce.g = upper(sce.g);


if prepare_input_only
    try
        fw = gui.myWaitbar(FigureHandle);
        run.py_panhumanpy(sce, wkdir, true, prepare_input_only);
        gui.myWaitbar(FigureHandle, fw);
        if strcmp(gui.myQuestdlg(FigureHandle, 'Input files prepared. Open the working folder?'),'Yes')
            winopen(wkdir);
        end
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    needupdatesce = false;
else
    fw = gui.myWaitbar(FigureHandle);
    try
        [c, T] = run.py_panhumanpy(sce, wkdir, true);
        assert(sce.NumCells==numel(c));
        stashname = pkg.i_stashcelltypehistory(sce);
        sce.c_cell_type_tx = c;
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    gui.myGuidata(FigureHandle, sce, src);
    needupdatesce = true;
    gui.myWaitbar(FigureHandle, fw);

    % Say what was assigned and where the labels it replaced went; the
    % pre-flight warning could only name the attribute generically.
    msg = sprintf('%d cell type(s) assigned to %d cells by panhumanpy.', ...
        numel(unique(string(c))), sce.NumCells);
    gui.myHelpdlg(FigureHandle, msg + gui.i_stashnotice(stashname));
end
end

function i_restoregenelist(sce, g)
% Put back the gene list the callback upper-cased for panhumanpy.
sce.g = g;
end
