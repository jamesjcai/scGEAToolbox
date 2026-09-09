function [requirerefresh] = callback_DoubletDetection(src, ~)

requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~gui.gui_showrefinfo('Scrublet [PMID:30954476]', FigureHandle)
    return;
end

% Native MATLAB by default; the Python package is offered only when Python
% is already configured, since it is no longer needed to run this at all.
usepy = false;
if pkg.i_checkpython
    backend = gui.myQuestdlg(FigureHandle, 'Choose Scrublet backend:', '', ...
        {'MATLAB (native)', 'Python (scrublet)'}, 'MATLAB (native)');
    if isempty(backend), return; end
    usepy = strcmp(backend, 'Python (scrublet)');
end
if numel(unique(sce.c_batch_id)) > 1
    if ~strcmp(gui.myQuestdlg(FigureHandle, ...
            ['"When working with data from multiple ' ...
            'samples, run Scrublet on each sample ' ...
            'separately." Your data contains multiple ' ...
            'samples (cells with different c_batch_id). ' ...
            'Continue?'],''), 'Yes')
        return;
    end
end


if usepy
    extprogname = 'py_scrublet';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    if ~gui.i_setpyenv([], [], FigureHandle), return; end
end

methodtag = 'scrublet';
info = [];
fw = gui.myWaitbar(FigureHandle);
try
    if usepy
        [isDoublet, doubletscore] = run.py_scrublet_new(sce.X, wkdir);
    else
        [isDoublet, doubletscore, info] = sc_scrublet(sce.X);
    end
    if isempty(isDoublet) || isempty(doubletscore)
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, "Running Error.");
        return;
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

% What this run can honestly be reported as. SC_SCRUBLET ends in three
% different states that all used to arrive here looking alike, and the two
% it warns about it warns about on the command window, which an App
% Designer user never sees. PKG.E_DOUBLETCALLSUMMARY tells them apart:
% "nothing was tested" is not "no doublet detected", and a run in which
% most simulated doublets are indistinguishable from single cells is not a
% result to act on.
S = pkg.e_doubletcallsummary(logical(isDoublet(:)), info);

if S.IsWarning
    gui.myWarndlg(FigureHandle, S.Message, S.Title);
else
    gui.myHelpdlg(FigureHandle, S.Message, S.Title);
end

if ~S.OfferRemoval
    return;
end

if sce.NumCells == length(doubletscore)
    % The preview is a plain figure. WindowStyle='modal' would keep it above
    % every other window and block input to them, including the confirmation
    % below.
    tmpf_doubletdetection = figure();
    gui.i_stemscatter(sce.s, doubletscore);
    ax = gca;
    zlabel(ax, 'Doublet Score');
    title(ax, sprintf('Doublet Detection (%s)', methodtag))

    % Ask through the preview figure, not FigureHandle: anchored to the main
    % app the question is drawn inside that window and ends up underneath the
    % preview. A traditional figure parent gets a questdlg of its own instead,
    % which stacks above the plot and leaves it visible while deciding.
    answer = gui.myQuestdlg(tmpf_doubletdetection, ...
        sprintf("Remove %d doublets?", sum(isDoublet)));

    % Close the preview whatever the answer, so declining does not leave an
    % orphan figure behind, and so the message below is not covered by it.
    if pkg.i_isvalid(tmpf_doubletdetection)
        close(tmpf_doubletdetection);
    end

    if strcmp(answer, 'Yes')
        sce = sce.removecells(isDoublet);
        gui.myGuidata(FigureHandle, sce, src);
        gui.myHelpdlg(FigureHandle, 'Doublets deleted.');
        requirerefresh = true;
    end
end
end
