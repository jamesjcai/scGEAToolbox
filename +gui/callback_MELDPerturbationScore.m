function callback_MELDPerturbationScore(src, ~)

hFig = src.Parent.Parent;
if ~gui.gui_showrefinfo('MELD [PMID:33558698]'), return; end

[ok] = gui.i_confirmscript('Run MELD Perturbation Score (MELD)?', ...
    'py_MELD', 'python');
if ~ok, return; end


sce = guidata(hFig);
if numel(unique(sce.c_batch_id)) < 2
    warndlg('No batch effect (SCE.C_BATCH_ID is empty)');
    return;
end

if ~gui.i_setpyenv
    return;
end

%fw=gui.gui_waitbar;
try
    id = sce.c_batch_id;
    if ~isnumeric(id)
        id = grp2idx(sce.c_batch_id);
        id = id(:);
    end    
    [score, T] = run.py_MELD(sce.X, id);
    if isempty(score) || size(score, 1) ~= size(sce.X, 2)
        %gui.gui_waitbar(fw);
        errordlg("MELD Running Error");
        return;
    end
catch ME
    %gui.gui_waitbar(fw);
    errordlg(ME.message);
    rethrow(ME);
end
% gui.gui_waitbar(fw);

hFig = figure;
gui.i_gscatter3(sce.s, score(:, 2));
colorbar
% tb = findall(hFig, 'tag', 'FigureToolBar'); % get the figure's toolbar handle
tb = uitoolbar('Parent', hFig);
pkg.i_addbutton2fig(tb, 'on', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');
gui.gui_3dcamera(tb, 'MELD_Scores');


if ~(ismcc || isdeployed)
    labels = {'Save score values to variable named:', 'Save score table to variable named:'};
    vars = {'MELDScores', 'MELDTable'};
    values = {score, T};
    export2wsdlg(labels, vars, values);
else
    gui.i_exporttable(T, false, 'MELDTable');
end

end
