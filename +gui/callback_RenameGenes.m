function [requirerefresh] = callback_RenameGenes(src)
requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = questdlg('Select genes to be renamed?');
if ~strcmp(answer, 'Yes'), return; end
[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    helpdlg('No gene selected.', '');
    return;
end
answer = questdlg('Paste new gene names?');
if ~strcmp(answer, 'Yes'), return; end
renamedglist = gui.i_inputgenelist(glist);

if length(glist) ~= length(renamedglist)
    % helpdlg('____.','');
    return;
end

[y, idx] = ismember(upper(glist), upper(sce.g));
if ~all(y)
    errordlg('Unspecific running error.');
    return;
end
sce.g(idx) = renamedglist;
requirerefresh = true;
helpdlg(sprintf('Renamed %d genes.', length(glist)), '');

guidata(FigureHandle, sce);
end
