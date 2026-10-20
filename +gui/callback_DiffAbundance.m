function callback_DiffAbundance(src, ~)
%CALLBACK_DIFFABUNDANCE  Analyze > Differential Abundance (DA) Analysis.
%
%   Asks for the cell-type variable, the condition variable and the
%   conditions to compare, and whether cells come from several replicate
%   samples, then runs SC_DIFFABUNDANCE: propeller with replicates, a
%   cell-level test without. Read-only: the dataset is not changed, only
%   its Methods text is recorded.

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('propeller [PMID:36005887]', FigureHandle), return; end

[ctype, ctLabel] = gui.i_select1class(sce, false, ...
    'Select the cell type (or cluster) variable:', 'Cell Type', FigureHandle);
if isempty(ctype), return; end
[cond, condLabel] = gui.i_select1class(sce, false, ...
    'Select the variable that separates the conditions:', 'Batch ID', ...
    FigureHandle);
if isempty(cond), return; end
if isequal(string(condLabel), string(ctLabel))
    gui.myWarndlg(FigureHandle, ['The condition variable must differ ' ...
        'from the cell type variable.']);
    return;
end

[picked, levels] = gui.i_selectgroupsubset(cond, condLabel, FigureHandle);
if isempty(picked), return; end
if numel(levels) < 2
    gui.myWarndlg(FigureHandle, 'Select at least two conditions to compare.');
    return;
end

answer = gui.myQuestdlg(FigureHandle, ['Do the cells come from several ' ...
    'biological samples (replicates) per condition? A replicate-aware ' ...
    'test needs a variable naming each cell''s sample.'], '', ...
    {'Yes, select sample variable', 'No replicates'}, ...
    'Yes, select sample variable');
if isempty(answer), return; end

sid = [];
sidLabel = "";
if startsWith(answer, 'Yes')
    [sid, sidLabel] = gui.i_select1class(sce, false, ...
        'Select the variable that identifies each cell''s sample:', ...
        'Batch ID', FigureHandle);
    if isempty(sid), return; end
    sidLabel = string(sidLabel);
    if sidLabel == string(condLabel)
        gui.myWarndlg(FigureHandle, ['The sample variable must differ ' ...
            'from the condition variable: each sample is one replicate ' ...
            'within a condition.']);
        return;
    end
end

% Cells without a label in any chosen variable cannot be counted.
ctype = i_tostring(ctype);
cond = i_tostring(cond);
keep = picked(:) & ctype ~= "" & cond ~= "";
if ~isempty(sid)
    sid = i_tostring(sid);
    keep = keep & sid ~= "";
    sid = sid(keep);
end
ctype = ctype(keep);
cond = cond(keep);
groups = unique(cond)';     % the order sc_diffabundance tests them in
if numel(groups) < 2
    gui.myWarndlg(FigureHandle, 'The selected conditions have no labelled cells to compare.');
    return;
end

[methodid, transform] = i_askmethod(sid, cond, groups, FigureHandle);
if isempty(methodid), return; end

fw = gui.myWaitbar(FigureHandle);
% The rank test's "underpowered" warning is shown as a dialog, not
% printed, so it is silenced here and read back from LASTWARN.
[lastMsg, lastId] = lastwarn();
lastwarn('');
ws = warning('off', 'sc_diffabundance:underpowered');
try
    [T, Counts] = sc_diffabundance(ctype, cond, sid, methodid, [], transform);
    warning(ws);
    [warnMsg, warnId] = lastwarn();
    lastwarn(lastMsg, lastId);
catch ME
    warning(ws);
    lastwarn(lastMsg, lastId);
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

[ptext, prefs] = pkg.i_damethodstext(T, CellTypeVar=string(ctLabel), ...
    GroupVar=string(condLabel), SampleVar=sidLabel, ...
    Groups=groups, NumSamples=numel(unique(sid)));
pkg.i_methodsreport(sce, "da", ptext, prefs);

if warnId == "sc_diffabundance:underpowered"
    gui.myWarndlg(FigureHandle, warnMsg);
end

i_showresult(T, Counts, groups, string(condLabel), FigureHandle);
end


function s = i_tostring(x)
% Per-cell labels as a string column; numeric IDs become their digits.
s = string(x(:));
s(ismissing(s)) = "";
end


function [methodid, transform] = i_askmethod(sid, cond, groups, parentfig)
% Offer the tests the design supports, the recommended one first.
methodid = [];
transform = "logit";
if ~isempty(sid)
    perGroup = arrayfun(@(g) numel(unique(sid(cond == g))), groups);
    if any(perGroup < 2)
        gui.myWarndlg(parentfig, sprintf(['Condition "%s" has only one ' ...
            'sample. A replicate-aware test needs at least two samples ' ...
            'in every condition.'], groups(find(perGroup < 2, 1))));
        return;
    end
    items = ["propeller, logit transform (recommended)", ...
        "propeller, arcsine square-root transform (less sensitive to outlier samples)", ...
        "Rank test on sample proportions (Wilcoxon / Kruskal-Wallis)"];
    codes = [5, 5, 4];
    trans = ["logit", "asin", "logit"];
    prompt = sprintf(['%d samples. propeller tests the per-sample ' ...
        'proportions and has power with 2-3 samples per condition; the ' ...
        'rank test needs about 5 or more.'], numel(unique(sid)));
else
    items = "Cell-level permutation test";
    codes = 2;
    trans = "logit";
    if numel(groups) == 2
        items(end+1) = "Fisher's exact test";
        codes(end+1) = 3;
        trans(end+1) = "logit";
    end
    prompt = ['Without replicates every cell is treated as independent, ' ...
        'so p-values are overconfident when samples differ.'];
end
[indx, tf] = gui.myListdlg(parentfig, cellstr(items), prompt, 1, false);
if ~tf || isempty(indx), return; end
methodid = codes(indx);
transform = trans(indx);
end


function i_showresult(T, Counts, groupNames, condLabel, parentfig)
% Grouped bars of each cell type's share per condition, significant cell
% types starred, with the full table one click away.
propVars = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, "Prop_"));
P = T{:, propVars};

hx = gui.myFigure(parentfig);
ax = hx.AxHandle;
b = bar(ax, categorical(T.CellType, T.CellType), P, 'grouped');
ylabel(ax, 'Proportion of cells');
legend(ax, groupNames, 'Interpreter', 'none', 'Location', 'best');
title(ax, sprintf('Cell type proportions by %s', condLabel), 'Interpreter', 'none');
subtitle(ax, extractAfter(string(T.Properties.Description), ": ") + ...
    "   * FDR < 0.05", 'Interpreter', 'none');
sig = T.FDR < 0.05;
if any(sig)
    top = max(P, [], 2);
    text(ax, find(sig), top(sig), '*', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 14);
end
ax.YLim(2) = ax.YLim(2)*1.08;
for k = 1:numel(b)
    b(k).BarWidth = 0.9;
end

hx.addCustomButton('off', @in_showtable, "table.gif", 'Show result table');
hx.addCustomButton('off', @in_savetable, "floppy-disk-arrow-in.jpg", ...
    'Save result table');
hx.show(parentfig);

    function in_showtable(~, ~)
        gui.i_viewtable({T, Counts}, hx.FigHandle, ["Test results", "Cell counts"]);
    end

    function in_savetable(~, ~)
        gui.i_exporttable(T, false, 'DATable', 'DiffAbundance', [], [], hx.FigHandle);
    end
end
