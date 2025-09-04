function callback_SelectCellsByClass(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

[ptsSelected] = gui.i_select1classcells(sce, true, FigureHandle);
if isempty(ptsSelected), return; end

parentax = findall(FigureHandle,'type','axes');

[ax, bx] = view(parentax);
fw = gui.myWaitbar(FigureHandle);
try
    scex = copy(sce).selectcells(ptsSelected); %#OK
    % scex.c=cLi(ci(idx));
    scex.c = sce.c(ptsSelected);
    if isa(src, 'matlab.apps.AppBase')
        a = scgeatoolApp(scex);
        view(a.UIAxes, [ax, bx]);
    else
        scgeatool(scex);
        view(ax, bx);
    end

    
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
end
