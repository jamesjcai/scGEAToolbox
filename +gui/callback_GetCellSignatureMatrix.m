function callback_GetCellSignatureMatrix(src, ~)
%CALLBACK_GETCELLSIGNATUREMATRIX Score cells for several signatures; radar plot.
%
%   The dialogs come in the order the analysis is framed: which cells to
%   work on (all, or chosen cell types/groups), which signatures, which
%   scoring method, and which groups to compare. Everything is asked before
%   anything is computed, so the one long wait comes at the end. The score
%   matrix is saved from the radar plot's toolbar rather than offered here.

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);

% SELECTCELLS mutates a handle object in place, so work on a copy or the
% subsetting below would remove cells from the app's own data.
sce = copy(sce_ori);

% ---- 1. Which cells to work on -------------------------------------------
answer = gui.myQuestdlg(FigureHandle, ['Score all cells, or only ' ...
    'selected cell types/groups?'], 'Cells to Analyze', ...
    {'All Cells', 'Select Cells', 'Cancel'}, 'All Cells');
switch answer
    case 'All Cells'
    case 'Select Cells'
        [thisc0, clabel0] = gui.i_selectnclass(sce, false, [], [], FigureHandle);
        if isempty(thisc0), return; end
        picked = gui.i_selectgroupsubset(thisc0, clabel0, FigureHandle);
        if isempty(picked) || ~any(picked), return; end
        sce.selectcells(picked);
    otherwise
        return;
end

% ---- 2. Which signatures --------------------------------------------------
[~, T] = pkg.e_cellscores([], [], 0);

% Same collection chooser as gui.callback_CompareCellScoreBtwCls - see
% gui.i_picksignaturetags.
[rows, taglabel] = gui.i_picksignaturetags(T, FigureHandle);
if isempty(rows), return; end

listitems = natsort(T.ScoreType(rows));
if strlength(taglabel) > 0
    % A named collection is the set the user asked for, so start with all
    % of it selected; they can still deselect.
    scoreprompt = char("Select Scores (" + taglabel + ")");
    preselected = true(numel(listitems), 1);
else
    scoreprompt = 'Select Scores';
    preselected = false(numel(listitems), 1);
end

if gui.i_isuifig(FigureHandle)
    [indx2, tf2] = gui.myListdlg(FigureHandle, listitems, ...
        scoreprompt, listitems(preselected), true);
else
    initval = find(preselected);
    if isempty(initval), initval = 1; end   % listdlg wants a valid index
    [indx2, tf2] = listdlg('PromptString', scoreprompt, ...
        'SelectionMode', 'multiple', 'ListString', ...
        listitems, 'ListSize', [260, 300], ...
        'InitialValue', initval);
end
if tf2 ~= 1 || isempty(indx2), return; end
scorenames = listitems(indx2);

[~, methodid] = gui.i_pickscoremethod([], FigureHandle);
if isempty(methodid), return; end

% ---- 3. Which groups to compare --------------------------------------------
answer = gui.myQuestdlg(FigureHandle, ...
    ['Compare signature scores between cell groups?', newline, newline, ...
    'Yes: one radar polygon per cell group.', newline, ...
    'No: a single radar polygon for all cells combined.'], ...
    'Cell group comparison');
switch answer
    case 'Yes'
        % Several grouping variables may be picked; they cross into one
        % composite label per cell ("Macrophages | IL").
        [thisc, clabel] = gui.i_selectnclass(sce, false, [], [], FigureHandle);
        if isempty(thisc), return; end
        if ~isscalar(unique(thisc))
            picked = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
            if isempty(picked) || ~any(picked), return; end
            sce.selectcells(picked);
            thisc = thisc(picked);
        end
    case 'No'
        thisc = repmat("All cells combined", sce.NumCells, 1);
    otherwise
        return;
end

% ---- 4. Score ---------------------------------------------------------------
n = numel(scorenames);
Y = zeros(sce.NumCells, n);
valid = false(n, 1);
fw = gui.myWaitbar(FigureHandle);
for k = 1:n
    a = scorenames{k};
    gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('Processing %s', a), (k - 1)/n);
    % A signature with too few expressed genes errors in e_cellscores; skip
    % it with a warning instead of aborting the whole analysis.
    try
        y = pkg.e_cellscores(sce.X, sce.g, a, methodid, false);
    catch ME
        warning('Skipping score "%s": %s', a, ME.message);
        y = [];
    end
    if ~isempty(y)
        Y(:, k) = y(:);
        valid(k) = true;
    end
end
gui.myWaitbar(FigureHandle, fw);

% Drop scores that could not be computed so one sparse signature does not
% abort the run or add an all-zero axis to the plots.
if ~any(valid)
    gui.myWarndlg(FigureHandle, ['No scores could be computed. The selected ' ...
        'signatures have too few expressed genes in this dataset.']);
    return;
end
if ~all(valid)
    gui.myWarndlg(FigureHandle, sprintf( ...
        '%d of %d score(s) skipped (too few expressed genes): %s', ...
        sum(~valid), n, strjoin(string(scorenames(~valid)), ', ')));
    Y = Y(:, valid);
    scorenames = scorenames(valid);
    n = numel(scorenames);
end

% ---- 5. Plot ----------------------------------------------------------------
labelx = scorenames';
if n >= 3
    % The radar plot's toolbar saves the score matrix.
    gui.i_spiderplot(Y, thisc, labelx, sce, FigureHandle);
    return;
end

% One or two scores make no radar, and a violin plot has no save button,
% so offer the matrix here instead.
c_cellid = string(sce.c_cell_id);
Tout = array2table(Y, 'VariableNames', scorenames, 'RowNames', ...
    matlab.lang.makeUniqueStrings(c_cellid));
Tout.Properties.DimensionNames{1} = 'Cell_ID';
gui.i_exporttable(Tout, true, 'Tcellsignmt', 'CellSignatTable', ...
    [], [], FigureHandle);
% One window, a tab per score, rather than a window per score.
labelinfo = struct('ylabel', 'Cellular score', 'xlabel', 'Cell group');
gui.sc_uitabgrpfig_vioplot(num2cell(Y, 1), string(labelx), thisc, ...
    FigureHandle, [], labelinfo);
end
