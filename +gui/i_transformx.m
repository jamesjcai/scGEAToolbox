function [X] = i_transformx(X, donorm, methodid, parentfig)

if nargin < 4, parentfig = []; end
if nargin < 3 || isempty(methodid), methodid = 3; end
if nargin < 2 || isempty(donorm), donorm = false; end
if nargin < 1
    X = nbinrnd(20, 0.98, 1000, 200);
    disp('Using simulated X.');
end
if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end
if donorm
    defaultans = 'Yes';
else
    defaultans = 'No';
end
answer = gui.myQuestdlg(parentfig, 'Normalize, transform or impute X? Select No to use untouched X', ...
'', {'Yes', 'No', 'Cancel'}, defaultans);
if strcmp(answer, 'Yes')

elseif strcmp(answer, 'No')
    return;
elseif strcmp(answer, 'Cancel')
    X = [];
    return;
else
    error('Wrong option');
end

listitems = {'(a): Library Size Normalization', ...
'(b): Log(x+1) Transformation', ...
'(c): (a) and (b)', ...
'(d): Shifted CLR PFlog1pPF Normalization', ...
'(e): DeSeq Normalization', ...
'(f): Pearson Residuals Transformation', ...
'(g): kNN Smoothing Transformation', ...
'(h): Freeman-Tukey Transformation', ...
'(i): MAGIC Imputation', ...
'(j): SCTransform (R/Seurat, vst v2)', ...
'(k): SCTransform (native MATLAB)'};
if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, listitems, ...
        'Select Method', listitems(methodid));
else
    [indx, tf] = listdlg('PromptString', {'Select Method'}, ...
        'SelectionMode', 'single', ...
        'ListString', listitems, 'ListSize', [220, 300], ...
        'InitialValue', methodid);
end

if tf == 1
    fw = gui.myWaitbar(parentfig);
    try
        switch indx
            case 1
                X = sc_norm(X, 'type', 'libsize');
            case 2
                X = log1p(X);
            case 3
                X = sc_norm(X, 'type', 'libsize');
                X = log1p(X);
            case 4
                X = sc_norm(X, 'type', 'shiftedclr');
            case 5
                X = sc_norm(X, 'type', 'deseq');
            case 6
                X = sc_transform(X, 'type', 'PearsonResiduals');
            case 7
                X = sc_transform(X, 'type', 'kNNSmoothing');
            case 8
                X = sc_transform(X, 'type', 'FreemanTukey');
            case 9
                X = run.ml_MAGIC(X, true);
            case 10
                X = sc_transform(X, 'type', 'SCTransform');
            case 11
                X = sc_transform(X, 'type', 'SCTransformMATLAB');
        end
    catch ME
        gui.myWaitbar(parentfig, fw);
        gui.myErrordlg(parentfig, ME.message)
        rethrow(ME)
    end
    gui.myWaitbar(parentfig, fw);
end

end
