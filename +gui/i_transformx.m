function [X] = i_transformx(X, donorm)

if nargin < 2, donorm = false; end

if nargin < 1
    X = nbinrnd(20, 0.98, 1000, 200);
    disp('Using simulated X.');
end
if donorm
    answer = questdlg('Normalize, transform or impute X?', ...
        '', 'Yes', 'No', 'Cancel', 'Yes');
else
    answer = questdlg('Normalize, transform or impute X?', ...
        '', 'Yes', 'No', 'Cancel', 'No');
end

if strcmp(answer, 'Yes')
    %
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
[indx, tf] = listdlg('PromptString', {'Select Method'}, ...
    'SelectionMode', 'single', ...
    'ListString', listitems, 'ListSize', [200, 300], ...
    'InitialValue', 3);
if tf == 1
    fw = gui.gui_waitbar;
    try
        switch indx
            case 1
                X = sc_norm(X);
            case 2
                X = log(X+1);
            case 3
                X = sc_norm(X);
                X = log(X+1);
            case 4
                X = sc_norm(X, 'type', 'deseq');
            case 5
                X = sc_transform(X, 'type', 'PearsonResiduals');
            case 6
                X = sc_transform(X, 'type', 'kNNSmoothing');
            case 7
                X = sc_transform(X, 'type', 'FreemanTukey');
            case 8
                X = run.mt_MAGIC(X, true);
        end
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message)
        rethrow(ME)
    end
    gui.gui_waitbar(fw);
end
end
