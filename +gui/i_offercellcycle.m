function isdone = i_offercellcycle(app)
%I_OFFERCELLCYCLE  Offer to estimate cell cycle phases the app shows as undetermined.
%
%   isdone = gui.i_offercellcycle(app)
%
% Called by Label Cell Groups after it has labelled the cells by Cell
% Cycle Phase. When every cell reads "undetermined" the phases were never
% estimated, and the plot - already showing that one group - is the
% reason to ask. Yes estimates the phases, keeps the old ones for Undo,
% and relabels the plot by the new phases. No leaves the plot as it is.
%
% ISDONE is true when the phases were estimated and the plot relabelled.
%
% See also GUI.CALLBACK_CELLCYCLEPOTENCY, GUI.I_SELECT1CLASS.

isdone = false;
sce = app.sce;
if ~isempty(sce.c_cell_cycle_tx) && ~all(strcmpi(string(sce.c_cell_cycle_tx), "undetermined"))
    return;
end

% Let the undetermined labels reach the screen, and stay there long enough
% to be read, before the dialog covers them.
drawnow;
pause(0.8);
answer = gui.myQuestdlg(app.UIFigure, ['Cell cycle phase is undetermined ' ...
    'for every cell. Estimate Cell Cycle Phase now?']);
if ~strcmp(answer, 'Yes'), return; end

gui.i_snapshot(app, 'Estimate Cell Cycle Phase');
fw = gui.myWaitbar(app.UIFigure);
try
    sce.estimatecellcycle(true, 1);
catch ME
    gui.myWaitbar(app.UIFigure, fw, true);
    gui.myErrordlg(app.UIFigure, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(app.UIFigure, fw);

% Relabel by the new phases, as Estimate Cell Cycle Phase does.
[app.c, app.cL] = findgroups(string(sce.c_cell_cycle_tx));
app.sce.c = app.c;
app.in_RefreshAll(true, false);
app.ix_labelclusters(true);
isdone = true;
end
