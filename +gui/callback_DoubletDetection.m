function [requirerefresh] = callback_DoubletDetection(src, ~)

requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~pkg.i_checkpython
    gui.myWarndlg(FigureHandle, 'Python not installed.');
    return;
end
if ~gui.gui_showrefinfo('Scrublet [PMID:30954476]', FigureHandle)
    return;
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


extprogname = 'py_scrublet';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
if ~gui.i_setpyenv([],[],FigureHandle), return; end


methodtag = 'scrublet';
try
    [isDoublet, doubletscore] = run.py_scrublet_new(sce.X, wkdir);
    if isempty(isDoublet) || isempty(doubletscore)
        gui.myErrordlg(FigureHandle, "Running Error.");
        return;
    end
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

if ~any(isDoublet)
    gui.myHelpdlg(FigureHandle, 'No doublet detected.');
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
