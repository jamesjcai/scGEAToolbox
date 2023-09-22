function xcallback_GetSTRINGNetworkImg(src, ~, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[speciestag] = i_selectspecies(3);