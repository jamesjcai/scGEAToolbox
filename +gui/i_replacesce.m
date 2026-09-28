function i_replacesce(app, sce, numCellsBefore)
%I_REPLACESCE Show SCE in the app's window in place of the current dataset.
%   gui.i_replacesce(app, sce, numCellsBefore) makes SCE the dataset of the
%   scgeatoolApp APP, redraws it keeping the view angle, and says how to get
%   the NUMCELLSBEFORE cells back.
%
%   This is where a subset goes instead of into a second scgeatoolApp
%   window: two main windows over each other, each with its own plots on
%   top, made it hard to tell which window a result belonged to. The
%   caller takes the Undo snapshot (gui.i_snapshot) before the data
%   changes, so Edit > Undo brings the full dataset back.
%
%   See also GUI.I_SNAPSHOT, GUI.CALLBACK_UNDO.

app.sce = sce;
% APP.C indexes the old cells; IN_REFRESHALL only recomputes it when empty.
[app.c, app.cL] = findgroups(string(sce.c));
app.in_RefreshAll(true, false);
gui.myHelpdlg(app.UIFigure, sprintf(['Now working on %d of the %d cells, ', ...
    'in this window.\n\nEdit > Undo (Ctrl+Z) brings back all %d.'], ...
    sce.NumCells, numCellsBefore, numCellsBefore));
end
