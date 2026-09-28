function callback_Dotplot(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL"). Downstream treats thisc as a
% per-cell label vector, so the composite needs no special handling.
[thisc, ~] = gui.i_selectnclass(sce,[],[],[],FigureHandle);
if isempty(thisc), return; end

[c, cL, noanswer] = gui.i_reordergroups(thisc, [], FigureHandle);
if noanswer, return; end

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
    return;
end
[Xt] = gui.i_transformx(sce.X, true, "libsize_log1p", FigureHandle);
if isempty(Xt), return; end
glist = glist(end:-1:1);

try
        gui.i_dotplot(Xt, sce.g, c, cL, glist, true, 'Dotplot', FigureHandle);
    catch ME
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    end
end
