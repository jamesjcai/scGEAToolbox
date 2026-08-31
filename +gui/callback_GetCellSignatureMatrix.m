function callback_GetCellSignatureMatrix(src, ~)
%     answer = gui.myQuestdlg(FigureHandle, ['This function ' ...
%         'calculates selected signature scores for each ' ...
%         'cell. You will get a signature matrix for cells.' ...
%         ' Continue?'],'');
%     if ~strcmp(answer,'Yes'), return; end


[FigureHandle, sce] = gui.gui_getfigsce(src);

[~, T] = pkg.e_cellscores([], [], 0);

% Same collection chooser as gui.callback_CompareCellScoreBtwCls - see
% gui.i_picksignaturetags. It replaces a hand-rolled list that dropped
% untagged signatures, carried a "-----" separator row that silently
% aborted the callback when clicked, and only preselected a collection
% inside the full 145-item list rather than narrowing to it.
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
        scoreprompt, listitems(preselected));
else
    initval = find(preselected);
    if isempty(initval), initval = 1; end   % listdlg wants a valid index
    [indx2, tf2] = listdlg('PromptString', scoreprompt, ...
        'SelectionMode', 'multiple', 'ListString', ...
        listitems, 'ListSize', [260, 300], ...
        'InitialValue', initval);
end

if tf2 ~= 1 || isempty(indx2), return; end

n = length(indx2);
Y = zeros(sce.NumCells, n);
scorenames = listitems(indx2);      % names of the selected scores
valid = false(n, 1);

[~, methodid] = gui.i_pickscoremethod([], FigureHandle);
if isempty(methodid), return; end

fw = gui.myWaitbar(FigureHandle);
for k = 1:n
    a = scorenames{k};
    b = sprintf('Processing %s', a);
    if n ~= k
        gui.myWaitbar(FigureHandle, fw, false, '', b, k/n);
    else
        gui.myWaitbar(FigureHandle, fw, false, '', b, (k - 1)/n);
    end
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
c_cellid=sce.c_cell_id;
if ~isstring(c_cellid)
    c_cellid=string(c_cellid);
end
T = array2table(Y, 'VariableNames', ...
scorenames, 'RowNames', ...
matlab.lang.makeUniqueStrings(c_cellid));
T.Properties.DimensionNames{1} = 'Cell_ID';
needwait = true;
gui.i_exporttable(T, needwait,'Tcellsignmt','CellSignatTable', ...
[],[],FigureHandle);

        % "Tcellattrib","CellAttribTable"
        % "Tviolindata","ViolinPlotTable"
        % "Tcrosstabul","CrosstabulTable"
        % "Tcellsignmt","CellSignatTable"

% assignin('base','Y',Y);
% assignin('base','listitems',listitems(indx2));
% assignin('base','labelx',listitems(indx2));

labelx = scorenames';
% T=table(Y,'VariableNames', ...
%     matlab.lang.makeValidName(listitems(indx2)));

answer = gui.myQuestdlg(FigureHandle, ...
    ['Compare signature scores between cell groups?', newline, newline, ...
    'Yes: one radar polygon per cell group.', newline, ...
    'No: a single radar polygon for all cells combined.'], ...
    'Cell group comparison');
if strcmp(answer, 'No')
    % Do not compare groups: treat all cells as one group so the radar plot
    % (n >= 3 scores) still shows the overall signature profile.
    thisc = repmat("All cells combined", sce.NumCells, 1);
elseif strcmp(answer, 'Yes')
    allowunique = false;
    % Several grouping variables may be picked; they cross into one composite
    % label per cell ("Macrophages | IL"). Downstream treats thisc as a
    % per-cell label vector, so the composite needs no special handling.
    [thisc] = gui.i_selectnclass(sce, allowunique, [], [], FigureHandle);
    if isempty(thisc), return; end
else
    return;
end

if n == 1
        gui.i_violinplot(Y, thisc, labelx, true, [], [], FigureHandle);
        xlabel('Cell group');
        ylabel('Cellular score');
    elseif n == 2
        gui.i_violinplot(Y(:, 1), thisc, labelx{1}, true, [], [], FigureHandle);
        xlabel('Cell group');
        ylabel('Cellular score');
        gui.i_violinplot(Y(:, 2), thisc, labelx{2}, true, [], [], FigureHandle);
        xlabel('Cell group');
        ylabel('Cellular score');

    elseif n >= 3
        gui.i_spiderplot(Y, thisc, labelx, sce, FigureHandle);
    end
end
