function callback_MELDPerturbationScore(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('MELD [PMID:33558698]', FigureHandle), return; end

[ok] = gui.i_confirmscript('Run MELD Perturbation Score (MELD)?', ...
    'py_MELD', 'python');
if ~ok, return; end



if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, 'No batch effect (SCE.C_BATCH_ID is empty)');
    return;
end

if ~gui.i_setpyenv([],[],FigureHandle)
    return;
end

%fw=gui.myWaitbar(FigureHandle);
try
    id = sce.c_batch_id;
    if ~isnumeric(id)
        id = findgroups(sce.c_batch_id);
        id = id(:);
    end    
    [score, T] = run.py_MELD(sce.X, id);
    if isempty(score) || size(score, 1) ~= size(sce.X, 2)
        %gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, "MELD error");
        return;
    end
catch ME
    %gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    rethrow(ME);
end
% gui.myWaitbar(FigureHandle, fw);

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
