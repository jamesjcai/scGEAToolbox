function [X] = i_transformx(X, donorm, methodid, parentfig)
%I_TRANSFORMX Ask whether and how to normalise, transform or impute X.
%   X = gui.i_transformx(X, donorm, methodid, parentfig) asks whether to
%   transform X (DONORM sets the default answer) and, if so, which way,
%   with METHODID preselected. X is [] if the user cancels.
%
%   METHODID names an item by its Key in gui.i_transformitems, e.g.
%   "pearson_residuals" or "libsize_log1p" (the default). A position in the
%   list is still accepted, but a name keeps meaning the same item when
%   items are inserted; a position does not. Either is checked before any
%   dialog opens.
%
%   See also gui.i_transformitems.

if nargin < 4, parentfig = []; end
if nargin < 3 || isempty(methodid), methodid = "libsize_log1p"; end
if nargin < 2 || isempty(donorm), donorm = false; end
if nargin < 1
    X = nbinrnd(20, 0.98, 1000, 200);
    disp('Using simulated X.');
end

items = gui.i_transformitems();
methodid = i_resolveitem(items, methodid);
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
if donorm
    defaultans = 'Yes';
else
    defaultans = 'No';
end
answer = gui.myQuestdlg(parentfig, ['Normalize, transform or impute the ', ...
    'expression data first? Choose No to use the data unchanged.'], ...
    'Transform Data', {'Yes', 'No', 'Cancel'}, defaultans);
if strcmp(answer, 'Yes')

elseif strcmp(answer, 'No')
    return;
else
    % 'Cancel', or the dialog closed: myQuestdlg returns '' for that, which
    % used to reach error('Wrong option') uncaught.
    X = [];
    return;
end

listitems = cellstr(items.Label);
if gui.i_isuifig(parentfig)
    % ALLOWMULTI false. It defaults to true in GUI.MYLISTDLG, so this
    % branch let the user pick several methods while the LISTDLG branch
    % below passes 'SelectionMode', 'single'. The switch beneath can only
    % honour one: MATLAB compares a switch expression to each case with
    % ISEQUAL, so a two-element INDX matches nothing, no branch runs, and X
    % comes back untransformed with no indication that the selection was
    % ignored.
    [indx, tf] = gui.myListdlg(parentfig, listitems, ...
        'Transform Data', listitems(methodid), false, true, [300, 450], ...
        'Select how to normalize, transform or impute the expression data.');
else
    [indx, tf] = listdlg('PromptString', {'Select how to transform the expression data'}, ...
        'SelectionMode', 'single', ...
        'ListString', listitems, 'ListSize', [220, 300], ...
        'InitialValue', methodid);
end

if tf == 1
    fw = gui.myWaitbar(parentfig);
    try
        % By key, so the cases follow the items wherever they sit in the list.
        if ~isscalar(indx)
            gui.myWaitbar(parentfig, fw);
            error('gui:i_transformx:unknownMethod', ...
                'Method index %s does not name one of the %d transforms.', ...
                mat2str(indx), numel(listitems));
        end
        switch items.Key(indx)
            case "libsize"
                X = sc_norm(X, 'type', 'libsize');
            case "log1p"
                X = log1p(X);
            case "libsize_log1p"
                X = sc_norm(X, 'type', 'libsize');
                X = log1p(X);
            case "shiftedclr"
                X = sc_norm(X, 'type', 'shiftedclr');
            case "deseq"
                X = sc_norm(X, 'type', 'deseq');
            case "pearson_residuals"
                X = sc_transform(X, 'type', 'PearsonResiduals');
            case "knn_smoothing"
                X = sc_transform(X, 'type', 'kNNSmoothing');
            case "freeman_tukey"
                X = sc_transform(X, 'type', 'FreemanTukey');
            case "magic"
                X = run.ml_MAGIC(X, true);
            case "sctransform_r"
                X = sc_transform(X, 'type', 'SCTransform');
            case "sctransform_matlab"
                X = sc_transform(X, 'type', 'SCTransformMATLAB');
            otherwise
                gui.myWaitbar(parentfig, fw);
                error('gui:i_transformx:unknownMethod', ...
                    'Transform "%s" has no implementation here.', items.Key(indx));
        end
    catch ME
        gui.myWaitbar(parentfig, fw, true);
        gui.myErrordlg(parentfig, ME.message, ME.identifier);
        % Reported; now stop the caller the way a cancel does. The rethrow
        % this replaces went uncaught in every caller, so each failed
        % transform also dumped "Error while evaluating Menu Callback".
        X = [];
        return;
    end
    gui.myWaitbar(parentfig, fw);
else
    % Cancelled. X = [] is how this function signals that, and it is what
    % every caller tests: callback_CompareGeneNetwork,
    % callback_BuildGeneNetwork and callback_Dotplot all do
    % `if isempty(Xt), return; end`. There was no else here at all, so a
    % cancelled method dialog returned X exactly as passed in -- raw
    % untransformed counts -- and the analysis ran on them as though the
    % user had chosen a transform. The first dialog in this function
    % already sets X = [] on Cancel; this is the same contract.
    X = [];
end

end

function idx = i_resolveitem(items, methodid)
% METHODID as a position in ITEMS, from a Key or a position.
if isnumeric(methodid)
    idx = methodid;
    if ~(isscalar(idx) && idx == round(idx) && idx >= 1 && idx <= height(items))
        error('gui:i_transformx:unknownMethod', ...
            'Method index %s does not name one of the %d transforms.', ...
            mat2str(idx), height(items));
    end
else
    idx = find(items.Key == lower(string(methodid)));
    if isempty(idx)
        error('gui:i_transformx:unknownMethod', ...
            'Unknown transform "%s". Use one of: %s.', string(methodid), ...
            strjoin(items.Key, ', '));
    end
end
end
