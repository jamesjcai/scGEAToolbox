function [X] = i_transformx(X, donorm, methodid, parentfig)

if nargin < 4, parentfig = []; end
if nargin < 3 || isempty(methodid), methodid = 3; end
if nargin < 2 || isempty(donorm), donorm = false; end

if nargin < 1
    X = nbinrnd(20, 0.98, 1000, 200);
    disp('Using simulated X.');
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
    '(d): DeSeq Normalization', ...
    '(e): Pearson Residuals Transformation', ...
    '(f): kNN Smoothing Transformation', ...
    '(g): Freeman-Tukey Transformation', ...
    '(h): MAGIC Imputation'};
    if gui.i_isuifig(parentfig)
        [indx, tf] = gui.myListdlg(parentfig, listitems, 'Select Method');
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
                X = sc_norm(X);
            case 2
                X = log1p(X);
            case 3
                X = sc_norm(X);
                X = log1p(X);
            case 4
                X = sc_norm(X, 'type', 'deseq');
            case 5
                X = sc_transform(X, 'type', 'PearsonResiduals');
            case 6
                X = sc_transform(X, 'type', 'kNNSmoothing');
            case 7
                X = sc_transform(X, 'type', 'FreemanTukey');
            case 8
                X = run.ml_MAGIC(X, true);
        end
    catch ME
        gui.myWaitbar(parentfig, fw);
        gui.myErrordlg(parentfig, ME.message)
        rethrow(ME)
    end
    gui.myWaitbar(parentfig, fw);
end
end
