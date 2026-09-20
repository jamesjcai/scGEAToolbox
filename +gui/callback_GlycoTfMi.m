function callback_GlycoTfMi(src, ~)
%GUI.CALLBACK_GLYCOTFMI  Menu callback: mutual information between
%transcription factors and glycogenes.
%
%   Wraps GLY.TFMI, which scores every TF against every glycogene by the
%   mutual information of their single-cell expression. MI is used rather
%   than correlation because TF-target relations have thresholds and
%   saturation in them, and a monotone measure misses those.
%
% See also GLY.TFMI, GLY.ENZONTO, PKG.E_GETTFLIST, SC_TFACTIVITY.

[FigureHandle, sce] = gui.gui_getfigsce(src);

species = gui.i_selectspecies(2, true, FigureHandle);
if isempty(species), return; end

% MI is linear in cells and quadratic in the bin grid, so the subsample cap
% is the runtime knob. The default is ample for the default bin count; a
% user who raises it should know what it costs.
answer = gui.i_inputdlg(['Maximum cells to use (subsampled; MI is linear ', ...
    'in cells, and 3000 is ample for the default 5 bins):'], '3000', FigureHandle);
if isempty(answer), return; end
maxCells = str2double(answer{1});
if ~isfinite(maxCells) || maxCells < 50 || maxCells ~= fix(maxCells)
    gui.myErrordlg(FigureHandle, ...
        'Enter a whole number of 50 or more.', 'gly.tfmi');
    return;
end

fw = gui.myWaitbar(FigureHandle);
try
    [~, info] = gly.tfmi(sce.X, sce.g, Species=string(species), ...
        MaxCells=maxCells);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'gly.tfmi');
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(info.tfs) || isempty(info.glycogenes)
    gui.myWarndlg(FigureHandle, ['Not enough TFs or glycogenes survived the ', ...
        'detection filter to estimate mutual information on this dataset.']);
    return;
end

msg = [sprintf('%d TFs against %d glycogenes, on %d cells.', ...
        numel(info.tfs), numel(info.glycogenes), info.nCells)
       ""
       "MI ranks pairs; it does not orient them. A high-scoring pair says " + ...
       "the two move together in a way a correlation would miss, not that " + ...
       "the TF drives the glycogene."];
if isfield(info, 'auprc') && istable(info.auprc) && ~isempty(info.auprc)
    msg(end+1) = "";
    msg(end+1) = "The estimator validation table (AUPRC against known edges) " + ...
        "opens alongside the results - read it before trusting the ranking.";
end
gui.myHelpdlg(FigureHandle, msg, 'TF-glycogene mutual information');

if isfield(info, 'auprc') && istable(info.auprc) && ~isempty(info.auprc)
    gui.i_viewtable(info.auprc, FigureHandle);
end
gui.i_viewtable(info.topTFs, FigureHandle);
gui.i_exporttable(info.topPairs, true, 'Tglycotfmi', 'GlycoTfMiTable', ...
    [], [], FigureHandle);
end
