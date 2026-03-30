function callback_TFActivity(src, ~)
% CALLBACK_TFACTIVITY  Estimate TF activity per cell using DoRothEA regulons.
%
% Dedicated entry point for TF activity analysis. Bypasses the generic
% gene-set-collection picker in callback_CompareCellScoreBtwCls and goes
% directly to species / TF / method selection.
%
% Workflow:
%   1. Ask: violin comparison across groups, or scatter/heatmap?
%   2. (If comparison) Select cell grouping variable.
%   3. Select species (human / mouse).
%   4. Select TF(s) to analyse (from DoRothEA list).
%   5. Select scoring method via i_picktfactivitymethod.
%   6. Run sc_tfactivity → activity matrix (nTF × nCells).
%   7. Extract selected TF(s) and display results.

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
if isempty(sce_ori), return; end
sce = copy(sce_ori);

% ---- Step 1: display mode --------------------------------------------------
aa = 'Yes, compare scores (violinplot)';
bb = 'No, just show values (scatter/heatmap)';
answer2 = gui.myQuestdlg(FigureHandle, ...
    ['Compare TF activity scores between cell groups (violinplot), ' ...
     'or display values for all cells?'], '', {aa, bb}, aa);
switch answer2
    case aa,  showcomparision = true;
    case bb,  showcomparision = false;
    otherwise, return;
end

% ---- Step 2: cell grouping (comparison mode only) -------------------------
if showcomparision
    [thisc] = gui.i_select1class(sce, false, [], [], FigureHandle);
    if isempty(thisc), return; end
    if isscalar(unique(thisc))
        gui.myWarndlg(FigureHandle, 'All cells are in the same group. Cannot compare.');
        return;
    end
else
    thisc = ones(sce.NumCells, 1);
end

% ---- Step 3: species -------------------------------------------------------
speciestag = gui.i_selectspecies(2, false, FigureHandle);
if isempty(speciestag), return; end

% ---- Step 4: select TF(s) BEFORE computation ------------------------------
% Load TF list from DoRothEA so the user commits to a selection before the
% (potentially slow) sc_tfactivity call runs.
[T] = pkg.e_gettflist(speciestag);
tflist_db = unique(T.tf);

if gui.i_isuifig(FigureHandle)
    [indx_sel, tfok] = gui.myListdlg(FigureHandle, tflist_db, ...
        'Select transcription factor(s) to analyse:');
else
    [indx_sel, tfok] = listdlg( ...
        'PromptString', 'Select transcription factor(s) to analyse:', ...
        'SelectionMode', 'multiple', ...
        'ListString', tflist_db, ...
        'ListSize', [260, 400]);
end
if ~tfok || isempty(indx_sel), return; end
selected_tfs_db = tflist_db(indx_sel);   % TF name(s) chosen by user

% ---- Step 5: scoring method ------------------------------------------------
[~, methodid] = gui.i_picktfactivitymethod(FigureHandle);
if isempty(methodid), return; end

% ---- Step 6: compute TF activity -------------------------------------------
% sc_tfactivity always computes all TFs; selection filters the output.
% Method 6 (NMF) has an internal waitbar; all others need one here.
if methodid ~= 6, fw = gui.myWaitbar(FigureHandle); end
try
    [cs, tflist] = sc_tfactivity(sce.X, sce.g, [], speciestag, methodid);
catch ME
    if methodid ~= 6, gui.myWaitbar(FigureHandle, fw, true); end
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
if methodid ~= 6, gui.myWaitbar(FigureHandle, fw); end

if isempty(cs) || isempty(tflist)
    gui.myWarndlg(FigureHandle, 'No TF activity scores could be computed.');
    return;
end

% ---- Step 7: extract selected TF(s) and display ---------------------------
% Match user-selected TF names against the computed tflist (which may differ
% from tflist_db if ULM filtered out low-target TFs).
[found, indx] = ismember(upper(string(selected_tfs_db)), upper(string(tflist)));
missing = selected_tfs_db(~found);
if ~isempty(missing)
    gui.myWarndlg(FigureHandle, sprintf( ...
        'TF(s) not found in activity results (too few targets?):\n%s', ...
        strjoin(missing, ', ')));
end
indx = indx(found);
if isempty(indx), return; end

scores_sel   = cs(indx, :);          % nSelected × nCells
labels_sel   = tflist(indx);         % matched TF names

if showcomparision
    if isscalar(indx)
        gui.sc_uitabgrpfig_vioplot(scores_sel(1, :), ...
            char(labels_sel(1)), thisc, FigureHandle);
    else
        y_cell = num2cell(scores_sel, 2);
        gui.sc_uitabgrpfig_vioplot(y_cell, ...
            cellstr(labels_sel), thisc, FigureHandle);
    end
else
    if isscalar(indx)
        gui.i_stemscatterfig(sce, scores_sel(1, :)', [], ...
            char(labels_sel(1)), FigureHandle);
    else
        gui.i_scoreheatmap(scores_sel, cellstr(labels_sel), sce, FigureHandle);
    end
end
end
