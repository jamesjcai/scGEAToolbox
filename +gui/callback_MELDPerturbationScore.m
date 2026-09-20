function callback_MELDPerturbationScore(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('MELD [PMID:33558698]', FigureHandle), return; end

% MELD returns one likelihood column per condition and this figure shows a
% single column, so the comparison has to be stated rather than inferred.
% Until now the grouping was hard-wired to SCE.C_BATCH_ID and the column
% plotted was column 2 - whichever level the labels happened to sort to
% second. Neither the variable nor the direction was ever shown to the user.
[grp, ctrlName, condName, clabel] = i_askcomparison(sce, FigureHandle);
if isempty(grp), return; end

% Cells outside the comparison keep a "" label and take no part in the
% graph: a third condition left in would be smoothed into both densities.
keep = grp ~= "";
id = categorical(grp(keep), [ctrlName, condName]);
idnum = double(id);          % 1 = control, 2 = condition, by construction

% Native MATLAB by default; the Python package is offered only when Python
% is already configured, since it is no longer needed to run this at all.
usepy = false;
if pkg.i_checkpython
    backend = gui.myQuestdlg(FigureHandle, 'Choose MELD backend:', '', ...
        {'MATLAB (native)', 'Python (meld)'}, 'MATLAB (native)');
    if isempty(backend), return; end
    usepy = strcmp(backend, 'Python (meld)');
end

if usepy
    [ok] = gui.i_confirmscript('Run MELD Perturbation Score (MELD)?', ...
        'py_MELD', 'python', FigureHandle);
    if ~ok, return; end
    if ~gui.i_setpyenv([], [], FigureHandle)
        return;
    end
end

info = [];
fw = gui.myWaitbar(FigureHandle);
try
    if usepy
        [score, T] = run.py_MELD(sce.X(:, keep), idnum);
    else
        [score, T, info] = sc_meld(sce.X(:, keep), id);
    end
    if isempty(score) || size(score, 1) ~= nnz(keep)
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, "MELD error");
        return;
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    rethrow(ME);
end
gui.myWaitbar(FigureHandle, fw);

% Smoothing moves a likelihood away from the null value even when nothing is
% happening, so the observed spread only means something next to the spread
% shuffled labels produce.
if ~isempty(info) && all(isfinite(info.nullSd))
    observedSd = std(score(:, 2));
    excluded = '';
    if ~all(keep)
        excluded = sprintf(['\n\n%d of %d cells are outside this ' ...
            'comparison and were left out; their score is NaN.'], ...
            nnz(~keep), numel(keep));
    end
    gui.myHelpdlg(FigureHandle, sprintf( ...
        ['Likelihood of %s against %s, grouped by %s.\n\nIt ranges %.2f ' ...
        'to %.2f. Its spread is %.3f, against %.3f from shuffled labels ' ...
        '(%.1fx). A cell is only notable when it sits well away from ' ...
        '%.2f, the share of the analysed cells %s contributes.%s'], ...
        condName, ctrlName, clabel, min(score(:, 2)), max(score(:, 2)), ...
        observedSd, info.nullSd(2), observedSd/max(info.nullSd(2), eps), ...
        mean(idnum == 2), condName, excluded), 'MELD');
end

hx = gui.myFigure(FigureHandle);
gui.i_gscatter3(sce.s(keep, :), score(:, 2), 1, 1, hx.AxHandle);
% A plain bar is correct now: I_GSCATTER3 passes a continuous vector to
% SCATTER as CData, so CLim is the likelihood's own range.
colorbar(hx.AxHandle);
title(hx.AxHandle, sprintf('MELD likelihood of %s (vs %s)', ...
    condName, ctrlName), 'Interpreter', 'none');
hx.show(FigureHandle);

% Exported values are padded back out to every cell, so they can be handed
% straight to SCE.SETCELLATTRIBUTE without the caller tracking which cells
% took part. Cells that did not are NaN.
scoreFull = nan(sce.NumCells, size(score, 2));
scoreFull(keep, :) = score;
TFull = array2table(scoreFull, "VariableNames", T.Properties.VariableNames);

if ~(ismcc || isdeployed)
    labels = {'Save score values to variable named:', 'Save score table to variable named:'};
    vars = {'MELDScores', 'MELDTable'};
    values = {scoreFull, TFull};
    export2wsdlg(labels, vars, values);
else
    gui.i_exporttable(TFull, false, 'MELDTable',[],[],[],hx.FigHandle);
end

end


function [grp, ctrlName, condName, clabel] = i_askcomparison(sce, parentfig)
%I_ASKCOMPARISON  Grouping variable, control group and condition group.
%
%   [grp, ctrlName, condName, clabel] = i_askcomparison(sce, parentfig)
%
%   GRP is a NumCells-by-1 string holding CTRLNAME or CONDNAME, and "" for
%   a cell that takes no part in the comparison - one carrying a third
%   condition's label, or no label at all. CLABEL names the variable the
%   labels came from, for the figure title.
%
%   An empty GRP means the user cancelled or no variable could supply two
%   groups. The dialogs are raised here, so the caller only has to check
%   for that.
%
%   MELD itself is symmetric: SC_MELD returns a likelihood per level, all
%   summing to 1 across a row, and needs no control. The control is asked
%   for because the figure shows one column. With two conditions the other
%   column is just 1 minus this one, so the choice is purely which way up
%   the perturbation reads.

grp = strings(0, 1);
ctrlName = "";
condName = "";

[thisc, clabel] = gui.i_select1class(sce, false, ...
    'Select the variable that separates the conditions:', 'Batch ID', ...
    parentfig);
if isempty(thisc), return; end
clabel = string(clabel);

lab = string(thisc(:));
lab(ismissing(lab)) = "";
lab = strip(lab);
levels = unique(lab(lab ~= ""));
if numel(levels) < 2
    gui.myWarndlg(parentfig, sprintf( ...
        ['MELD compares two conditions, but "%s" defines %d group. ' ...
        'Pick a variable with at least two levels.'], ...
        clabel, numel(levels)));
    return;
end

% PREFERSEL goes in as a char, not a 1-by-1 cellstr: MYLISTDLG hands it
% straight to UILISTBOX's Value, which wants a bare element of Items when
% MultiSelect is off.
[indx, tf] = gui.myListdlg(parentfig, cellstr(levels), 'Control group', ...
    char(levels(1)), false, true, [300, 300], sprintf( ...
    ['Which level of "%s" is the control (reference)? The score shown ' ...
    'is the likelihood of the other group, so this sets which way round ' ...
    'the perturbation reads.'], clabel));
if tf ~= 1, return; end
ctrlName = levels(indx);

rest = levels(levels ~= ctrlName);
if isscalar(rest)
    condName = rest;
    grp = strings(numel(lab), 1);
    grp(lab == ctrlName) = ctrlName;
    grp(lab == condName) = condName;
    return;
end

% The pooled entry is appended last and is identified by its index, not its
% text: a level genuinely named "All other groups (pooled)" would otherwise
% take the wrong branch.
items = [rest; "All other groups (pooled)"];
[indx, tf] = gui.myListdlg(parentfig, cellstr(items), 'Condition group', ...
    char(items(1)), false, true, [300, 300], sprintf( ...
    ['Which group is scored against "%s"? Picking one group restricts ' ...
    'the analysis to those two conditions; pooling keeps every labelled ' ...
    'cell and treats the rest as a single condition.'], ctrlName));
if tf ~= 1
    ctrlName = "";
    return;
end

grp = strings(numel(lab), 1);
grp(lab == ctrlName) = ctrlName;
if indx == numel(items)
    % "not WT" cannot collide with the control's own name, which a fixed
    % label such as "All others" could.
    condName = "not " + ctrlName;
    grp(lab ~= "" & lab ~= ctrlName) = condName;
else
    condName = items(indx);
    grp(lab == condName) = condName;
end

end
