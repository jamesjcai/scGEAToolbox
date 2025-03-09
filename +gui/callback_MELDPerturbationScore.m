function callback_MELDPerturbationScore(src, ~)

[FigureHandle, sce, isui] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('MELD [PMID:33558698]'), return; end

[ok] = gui.i_confirmscript('Run MELD Perturbation Score (MELD)?', ...
    'py_MELD', 'python');
if ~ok, return; end



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

hx=gui.myFigure;
gui.i_gscatter3(sce.s, score(:, 2));
colorbar
hx.show(FigureHandle);

if ~(ismcc || isdeployed)
    labels = {'Save score values to variable named:', 'Save score table to variable named:'};
    vars = {'MELDScores', 'MELDTable'};
    values = {score, T};
    export2wsdlg(labels, vars, values);
else
    gui.i_exporttable(T, false, 'MELDTable',[],[],[],FigureHandle);
end

end
