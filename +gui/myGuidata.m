function myGuidata(FigureHandle, sce, src)
    if isa(src, 'matlab.apps.AppBase')
        src.sce = sce;
    else
        guidata(FigureHandle, sce);
    end
end