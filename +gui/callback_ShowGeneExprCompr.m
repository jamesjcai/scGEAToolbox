function callback_ShowGeneExprCompr(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

[axx, bxx] = view(findall(FigureHandle,'type','axes'));

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist), return; end

allowunique = false;
% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL"). Downstream treats thisc as a
% per-cell label vector, so the composite needs no special handling.
[thisc] = gui.i_selectnclass(sce, allowunique,[],[],FigureHandle);
if isempty(thisc), return; end
% The comparison draws two groups side by side (sc_uitabgrpfig_expcomp
% uses c == 1 and c == 2), so exactly two are needed. One group used to
% leave SCE1 unset after a "Continue?" Yes; one or 3+ picked from the list
% crashed or silently dropped groups, both with the progress bar left open.
[ci, cLi] = findgroups(string(thisc));
if isscalar(cLi)
    gui.myWarndlg(FigureHandle, ['All cells are in the same group, so ' ...
        'there is nothing to compare. Choose a grouping with at least two groups.']);
    return;
end
listitems = natsort(cLi);
if gui.i_isuifig(FigureHandle)
    [indxx, tfx] = gui.myListdlg(FigureHandle, ...
        listitems, 'Select two groups:', ...
        listitems(1:2), true);
else
    [indxx, tfx] = listdlg('PromptString', ...
        {'Select two groups:'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listitems, ...
        'InitialValue', 1:2, 'ListSize', [220, 300]);
end
if tfx ~= 1, return; end
if numel(indxx) ~= 2
    gui.myWarndlg(FigureHandle, sprintf( ...
        'Select exactly two groups to compare; %d were selected.', numel(indxx)));
    return;
end
[y1, idx1] = ismember(listitems(indxx), cLi);
assert(all(y1));
idx2 = ismember(ci, idx1);
sce1 = copy(sce).selectcells(idx2);  % OK
thisc = thisc(idx2);

fw = gui.myWaitbar(FigureHandle);
try
    gui.sc_uitabgrpfig_expcomp(sce1, glist, FigureHandle, [axx, bxx], thisc);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

end
