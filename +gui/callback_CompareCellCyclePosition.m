function callback_CompareCellCyclePosition(src, ~)
%CALLBACK_COMPARECELLCYCLEPOSITION Compare cell cycle position between two groups.
%   Analyze > Estimate Cell States > Compare Cell Cycle Position Between
%   Groups... Picks two cell groups, the first taken as the reference
%   (control), and runs sc_cellcyclecompare on them: tricycle places every
%   cell on the cycle and permutation tests ask whether the second group
%   sits elsewhere (mean shift) or bunches up (concentration). Draws the
%   comparison figure and puts the result struct in the base workspace.
%
%   See also sc_cellcyclecompare, run.ml_tricycle.

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, ['Place the cells of two groups on the ' ...
    'cell cycle with tricycle and test whether the second group sits ' ...
    'elsewhere on it than the first, the reference. Continue?']);
if ~strcmp(answer, 'Yes'), return; end

[i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
if isscalar(i1) || isscalar(i2), return; end

speciestag = gui.i_selectspecies(2, false, FigureHandle);
if isempty(speciestag), return; end

groupNames = [string(cL1), string(cL2)];
% Classic figure: the plot uses polaraxes in a tiledlayout. It stays hidden
% until SHOW, so it can be drawn into while the work runs.
hx = gui.myFigure(FigureHandle, true);
delete(hx.AxHandle);
fw = gui.myWaitbar(FigureHandle);
try
    result = sc_cellcyclecompare(sce.X(:, i1), sce.g, sce.X(:, i2), sce.g, ...
        Species=speciestag, GroupNames=groupNames, Parent=hx.FigHandle);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    delete(hx.FigHandle);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
hx.show(FigureHandle);

assignin("base", "cellCycleCompareResult", result);
gui.myHelpdlg(FigureHandle, sprintf(['%s (n = %d) vs %s (n = %d): mean shift ' ...
    '%.2f pi (p = %.3g), concentration %.2f vs %.2f (p = %.3g). The result is ' ...
    'in the workspace as cellCycleCompareResult.'], groupNames(1), ...
    result.NumCells(1), groupNames(2), result.NumCells(2), result.MeanShift/pi, ...
    result.PValueMeanShift, result.Concentration(1), result.Concentration(2), ...
    result.PValueConcentration));
end
