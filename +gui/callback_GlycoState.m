function callback_GlycoState(src, ~)
%GUI.CALLBACK_GLYCOSTATE  Menu callback: per-cell glycobiological state.
%
%   Wraps GLY.STATE, which scores every cell against the curated
%   glycobiology modules (GLY.GENESETS). A cell's row of scores is its
%   glycobiological state - the relative activity of N-/O-glycosylation,
%   sialylation, fucosylation, glycosaminoglycan and glycosphingolipid
%   metabolism, glycan degradation and glycan recognition.
%
%   Scoring the same collection is also reachable through Analyze > Gene
%   Program (Cell Score) Analysis, where "Glycobiology" is one of the
%   gene-set collections. This entry exists because it scores the whole
%   collection in one step and reports how many of each module's genes were
%   actually found, which is what decides whether a module score means
%   anything on this dataset.
%
% See also GLY.STATE, GLY.GENESETS, GUI.CALLBACK_GLYCODETECT.

[FigureHandle, sce] = gui.gui_getfigsce(src);

[~, methodid] = gui.i_pickscoremethod([], FigureHandle);
if isempty(methodid), return; end

minGenes = in_askmingenes(FigureHandle);
if isempty(minGenes), return; end

fw = gui.myWaitbar(FigureHandle);
try
    [cs, setnames, ncommon] = gly.state(sce.X, sce.g, methodid, minGenes);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'gly.state');
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(setnames)
    gui.myWarndlg(FigureHandle, sprintf( ...
        ['No glyco module had at least %d of its genes present in this ', ...
         'dataset, so nothing could be scored. Lower the threshold, or ', ...
         'check that the gene symbols match the collection (human ', ...
         'symbols, uppercase).'], minGenes));
    return;
end

% How many of each module's genes were found is not decoration: a module
% scored off three genes out of forty is a different object from one scored
% off all forty, and the score alone does not say which it is.
gui.i_viewtable(table(setnames, ncommon, ...
    VariableNames = ["module", "genesFound"]), FigureHandle);

cellid = string(sce.c_cell_id);
% makeUniqueStrings around makeValidName: the curated collection's names
% do not collide today, but array2table errors outright if two ever reduce
% to the same identifier, and a rename in GLY.GENESETS should not be able
% to break this callback.
T = array2table(cs', RowNames = matlab.lang.makeUniqueStrings(cellid), ...
    VariableNames = matlab.lang.makeUniqueStrings( ...
        matlab.lang.makeValidName(setnames)));
T.Properties.DimensionNames{1} = 'Cell_ID';
gui.i_exporttable(T, true, 'Tglycostate', 'GlycoStateTable', ...
    [], [], FigureHandle);

% A heatmap, not a radar plot: the collection is two dozen modules, and a
% two-dozen-axis polygon is a shape rather than a reading. This is the same
% view GUI.CALLBACK_TFACTIVITY gives its score matrix, which is the closest
% existing analysis in the toolbox.
answer = gui.myQuestdlg(FigureHandle, [ ...
    "View the module scores as a heatmap?"
    ""
    "Cells are grouped by an attribute you pick, so each module reads as a " + ...
    "band across the groups."]);
if strcmp(answer, 'Yes')
    gui.i_scoreheatmap(cs, cellstr(setnames), sce, FigureHandle);
end
end


function minGenes = in_askmingenes(parentfig)
% A module whose genes are mostly absent still produces a score, and that
% score looks exactly like a real one. MINGENES is what keeps it out.

minGenes = [];
answer = gui.i_inputdlg(['Minimum genes of a module that must be present ', ...
    'in the data (modules below this are dropped):'], '3', parentfig);
if isempty(answer), return; end

v = str2double(answer{1});
if ~isfinite(v) || v < 1 || v ~= fix(v)
    gui.myErrordlg(parentfig, ...
        'Enter a whole number of 1 or more.', 'gly.state');
    return;
end
minGenes = v;
end
