function callback_SelectCellsByClass(src, ~)
% The selected cells replace the dataset in this window; the menu handler
% takes the Undo snapshot, so Ctrl+Z brings the rest back. They used to open
% in a second scgeatoolApp window, over this one.

[FigureHandle, sce] = gui.gui_getfigsce(src);

[ptsSelected] = gui.i_select1classcells(sce, true, FigureHandle);
if isempty(ptsSelected), return; end

parentax = findall(FigureHandle,'type','axes');

[ax, bx] = view(parentax);
fw = gui.myWaitbar(FigureHandle);
try
    scex = copy(sce).selectcells(ptsSelected); % OK
    scex.c = sce.c(ptsSelected);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isa(src, 'matlab.apps.AppBase')
    gui.i_replacesce(src, scex, sce.NumCells);
else
    scgeatool(scex);
    view(ax, bx);
end
end
