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
    tmpf_doubletdetection = figure('WindowStyle', 'modal');
    gui.i_stemscatter(sce.s, doubletscore);
    ax = gca;
    zlabel(ax, 'Doublet Score');
    title(ax, sprintf('Doublet Detection (%s)', methodtag))
    if strcmp(gui.myQuestdlg(FigureHandle, ...
            sprintf("Remove %d doublets?", sum(isDoublet))),'Yes')
        close(tmpf_doubletdetection);
        sce = sce.removecells(isDoublet);
        gui.myGuidata(FigureHandle, sce, src);
        gui.myHelpdlg(FigureHandle, 'Doublets deleted.');
        requirerefresh = true;
    end
end
end
