function [requirerefresh, s] = callback_MergeSCEs(src, ~)
% CALLBACK_MERGESCES - Merge the open dataset with other SCEs.
%
%   [requirerefresh, s] = gui.callback_MergeSCEs(src)
%
% Behind Edit > Merge Current Dataset with Others. The open dataset comes
% first; the others are SCE variables in the base workspace or SCE .mat
% files. The merged dataset replaces the open one, after a confirmation that
% says what will replace it.
%
% Merging datasets none of which is open yet is File > Import Data, which
% takes several SCE files or several workspace variables. Both go through
% GUI.I_MERGESCES, so they ask the same questions and label batches the same
% way. The second input is ignored; it once chose between the workspace and
% files, which is now asked here.
%
% When src is the app, the merged dataset is installed and the plot redrawn
% here, and REQUIREREFRESH comes back false: the menu handler in
% scgeatoolApp.mlapp would otherwise redraw a second time and follow up with
% a message built from the dataset names in upper case. For any other src,
% REQUIREREFRESH is true and S lists the merged datasets.
%
% see also: gui.i_mergesces, gui.i_pickworkspacesces, gui.i_loadscefiles

requirerefresh = false;
s = "";
[parentfig, cursce] = gui.gui_getfigsce(src);
if isempty(cursce) || cursce.NumCells == 0
    gui.myWarndlg(parentfig, ['There is no dataset open to merge with. Use ' ...
        'File > Import Data, which merges several SCE files or workspace ' ...
        'variables as it reads them.']);
    return;
end

answer = gui.myQuestdlg(parentfig, ['Merge the open dataset with SCE ' ...
    'variables in the base workspace, or with SCE data files?'], ...
    'Merge Datasets', {'Workspace Variables', 'SCE Data Files'}, ...
    'Workspace Variables');
switch answer
    case 'Workspace Variables'
        [others, names] = gui.i_pickworkspacesces(parentfig, ...
            'Select the datasets to merge with the open one.', cursce);
    case 'SCE Data Files'
        [fname, pathname] = uigetfile({'*.mat', 'SCE Data Files (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Select SCE Data Files to Merge with the Open Dataset', ...
            'MultiSelect', 'on');
        if isequal(fname, 0), return; end
        [others, names] = gui.i_loadscefiles(parentfig, pathname, fname);
    otherwise
        return;
end
if isempty(others), return; end

insce = [{cursce}, others];
names = ["current", names];
sce = gui.i_mergesces(parentfig, insce, names, @in_confirm);
if isempty(sce), return; end

gui.myGuidata(parentfig, sce, src);
if ~isa(src, 'matlab.apps.AppBase')
    requirerefresh = true;
    s = strjoin(names, ",");
    return;
end

src.sce = sce;
[src.c, src.cL] = findgroups(string(sce.c_batch_id));
src.sce.c = src.c;
src.in_RefreshAll(true, false);
gui.myHelpdlg(parentfig, sprintf(['Merged %d datasets into %d cells and ' ...
    '%d genes, colored by batch (%d batches).'], numel(insce), ...
    sce.NumCells, sce.NumGenes, numel(src.cL)));

    function tf = in_confirm(ncells)
        msg = sprintf(['The merged dataset (%d cells from %s) will replace ' ...
            'the one open in this window. Edit > Undo brings it back. ' ...
            'Continue?'], ncells, strjoin(names, ", "));
        tf = strcmp(gui.myQuestdlg(parentfig, msg, 'Merge Datasets', ...
            {'Merge', 'Cancel'}, 'Merge', 'warning'), 'Merge');
    end
end
