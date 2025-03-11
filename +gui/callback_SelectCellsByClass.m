function callback_SelectCellsByClass(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);


[ptsSelected] = gui.i_select1classcells(sce, true);
if isempty(ptsSelected), return; end

[ax, bx] = view(findall(FigureHandle,'type','axes'));
fw = gui.gui_waitbar;
try
    scex = selectcells(sce, ptsSelected);
    % scex.c=cLi(ci(idx));
    scex.c = sce.c(ptsSelected);
    scgeatool(scex);
    view(ax, bx);
catch ME
    gui.gui_waitbar(fw, true);
    errordlg(ME.message);
    return;
end
gui.gui_waitbar(fw);
end
