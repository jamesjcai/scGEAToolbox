function callback_SelectCellsByClass(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

[ptsSelected] = gui.i_select1classcells(sce, true, FigureHandle);
if isempty(ptsSelected), return; end

parentax = findall(FigureHandle,'type','axes')

[ax, bx] = view(parentax);
fw = gui.myWaitbar(FigureHandle);
try
    scex = selectcells(sce, ptsSelected);
    % scex.c=cLi(ci(idx));
    scex.c = sce.c(ptsSelected);
    a = scgeatoolApp(scex);
    view(a.UIAxes, [ax, bx]);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
end
