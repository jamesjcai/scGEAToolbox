function [requirerefresh] = callback_RenameGenes(src)
requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, 'Select genes to be renamed?');
if ~strcmp(answer, 'Yes'), return; end
[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
    return;
end
answer = gui.myQuestdlg(FigureHandle, 'Paste new gene names?');
if ~strcmp(answer, 'Yes'), return; end
renamedglist = gui.i_inputgenelist(glist, [], FigureHandle);

if isempty(renamedglist), return; end   % paste cancelled
if length(glist) ~= length(renamedglist)
    % Said nothing before, so a miscounted paste looked like a rename.
    gui.myWarndlg(FigureHandle, sprintf(['%s selected but %s pasted; ' ...
        'the two lists must be the same length. Nothing was renamed.'], ...
        pkg.i_plural(length(glist), 'gene'), ...
        pkg.i_plural(length(renamedglist), 'name')));
    return;
end

[y, idx] = ismember(upper(glist), upper(sce.g));
if ~all(y)
    gui.myErrordlg(FigureHandle, 'Unspecific running error.');
    return;
end
sce.g(idx) = renamedglist;
requirerefresh = true;
gui.myHelpdlg(FigureHandle, sprintf('Renamed %s.', ...
    pkg.i_plural(length(glist), 'gene')), '');

gui.myGuidata(FigureHandle, sce, src);
end
