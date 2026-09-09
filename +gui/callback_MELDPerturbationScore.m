function callback_MELDPerturbationScore(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('MELD [PMID:33558698]', FigureHandle), return; end

if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, 'No batch effect (SCE.C_BATCH_ID is empty)');
    return;
end

% Native MATLAB by default; the Python package is offered only when Python
% is already configured, since it is no longer needed to run this at all.
usepy = false;
if pkg.i_checkpython
    backend = gui.myQuestdlg(FigureHandle, 'Choose MELD backend:', '', ...
        {'MATLAB (native)', 'Python (meld)'}, 'MATLAB (native)');
    if isempty(backend), return; end
    usepy = strcmp(backend, 'Python (meld)');
end

if usepy
    [ok] = gui.i_confirmscript('Run MELD Perturbation Score (MELD)?', ...
        'py_MELD', 'python', FigureHandle);
    if ~ok, return; end
    if ~gui.i_setpyenv([], [], FigureHandle)
        return;
    end
end

info = [];
fw = gui.myWaitbar(FigureHandle);
try
    id = sce.c_batch_id;
    if ~isnumeric(id)
        id = findgroups(sce.c_batch_id);
        id = id(:);
    end
    if usepy
        [score, T] = run.py_MELD(sce.X, id);
    else
        [score, T, info] = sc_meld(sce.X, id);
    end
    if isempty(score) || size(score, 1) ~= size(sce.X, 2)
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, "MELD error");
        return;
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    rethrow(ME);
end
gui.myWaitbar(FigureHandle, fw);

% Smoothing moves a likelihood away from the null value even when nothing is
% happening, so the observed spread only means something next to the spread
% shuffled labels produce.
if ~isempty(info) && all(isfinite(info.nullSd))
    observedSd = std(score(:, 2));
    gui.myHelpdlg(FigureHandle, sprintf( ...
        ['Likelihood for %s ranges %.2f to %.2f.\n\nIts spread is %.3f, ' ...
        'against %.3f from shuffled labels (%.1fx). A cell is only ' ...
        'notable when it sits well away from %.2f, the share of cells ' ...
        'this sample contributes.'], ...
        string(info.levels(2)), min(score(:, 2)), max(score(:, 2)), ...
        observedSd, info.nullSd(2), observedSd/max(info.nullSd(2), eps), ...
        mean(id == 2)), 'MELD');
end

hx = gui.myFigure(FigureHandle);
gui.i_gscatter3(sce.s, score(:, 2), 1, 1, hx.AxHandle);
colorbar(hx.AxHandle);
hx.show(FigureHandle);

if ~(ismcc || isdeployed)
    labels = {'Save score values to variable named:', 'Save score table to variable named:'};
    vars = {'MELDScores', 'MELDTable'};
    values = {score, T};
    export2wsdlg(labels, vars, values);
else
    gui.i_exporttable(T, false, 'MELDTable',[],[],[],hx.FigHandle);
end

end
