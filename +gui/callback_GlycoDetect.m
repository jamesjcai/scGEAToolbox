function callback_GlycoDetect(src, ~)
%GUI.CALLBACK_GLYCODETECT  Menu callback: is this dataset deep enough for
%glyco analysis?
%
%   Wraps GLY.DETECT. Glycogenes are lowly expressed, so a glyco readout
%   sits further down the detection curve than a marker-gene analysis and
%   moves more for a given change in depth. Run this before GLY.STATE or
%   GLY.ENRICH: both rest on detection, and neither can tell you whether
%   the detection you have is the biology or the library.
%
% See also GLY.DETECT, GUI.CALLBACK_GLYCOENRICH, GUI.CALLBACK_GLYCOSTATE.

[FigureHandle, sce] = gui.gui_getfigsce(src);

% Detection against depth is meaningless once every cell has been scaled to
% the same library size. GLY.DETECT warns about this on the console, where
% a GUI user will not see it, so ask here instead.
if in_lookssizenormalized(sce.X)
    answer = gui.myQuestdlg(FigureHandle, [ ...
        "Every cell in this dataset has nearly the same library size, so " + ...
        "the counts look library-size normalized."
        ""
        "This analysis reads detection against sequencing depth, which " + ...
        "normalized input flattens - the curve and the saturation estimate " + ...
        "will not mean anything. Re-import the raw counts for a usable " + ...
        "answer."
        ""
        "Continue anyway?"]);
    if ~strcmp(answer, 'Yes'), return; end
end

% The per-group breakdown is the quickest way to see a depth imbalance
% that would confound GLY.ENRICH, so offer it, but do not insist: the
% ungrouped verdict is useful on its own.
grp = [];
answer = gui.myQuestdlg(FigureHandle, [ ...
    "Break the report down by cell group?"
    ""
    "Yes: also report median depth and glycogene detection per group, " + ...
    "and the ratio between the deepest and shallowest. That ratio is what " + ...
    "says whether a later enrichment comparison is confounded."
    ""
    "No: one verdict for the whole dataset."]);
switch answer
    case 'Yes'
        grp = gui.i_getcellgroups(sce, FigureHandle);
        if isempty(grp), return; end
    case 'No'
        % Ungrouped run - grp stays empty.
    otherwise
        return;
end

fw = gui.myWaitbar(FigureHandle);
try
    [T, info] = gly.detect(sce.X, sce.g, Group=grp, Plot=true);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'gly.detect');
    return;
end
gui.myWaitbar(FigureHandle, fw);

msg = [string(info.verdict)
       ""
       sprintf('Glycogenes present in the data: %d', info.nGlycoInData)
       sprintf('Median UMI per cell: %.0f', info.medianUMI)
       sprintf('Median genes detected per cell: %.0f', info.medianDetected)];
if isfinite(info.saturation)
    msg(end+1) = sprintf('Fraction of the reachable glycogene repertoire captured: %.0f%%', ...
        100 * info.saturation);
    msg(end+1) = sprintf('UMI needed for 90%% of the ceiling: %.0f', info.depthFor90);
end
if isfield(info, 'depthRatioDetected') && ~isempty(info.depthRatioDetected) ...
        && isfinite(info.depthRatioDetected)
    msg(end+1) = "";
    msg(end+1) = sprintf(['Across groups, the deepest is %.1fx the shallowest in ', ...
        'genes detected (%.1fx in UMI). It is the detected ratio that governs ', ...
        'whether a detection-based comparison is confounded.'], ...
        info.depthRatioDetected, info.depthRatioUMI);
end
gui.myHelpdlg(FigureHandle, msg, 'Glycogene detection depth');

if isfield(info, 'byGroup') && istable(info.byGroup) && ~isempty(info.byGroup)
    gui.i_viewtable(info.byGroup, FigureHandle);
end

gui.i_exporttable(T, true, 'Tglycodetect', 'GlycoDetectionTable', ...
    [], [], FigureHandle);
end


function tf = in_lookssizenormalized(X)
% The same relative test GLY.DETECT makes internally: library-size
% normalization leaves every column summing to the same 1e4-scale number up
% to float error, which is far too close together for an absolute epsilon.

s = full(double(sum(X, 1)));
tf = ~isempty(s) && (max(s) - min(s)) / max(max(s), eps) < 1e-6;
end
