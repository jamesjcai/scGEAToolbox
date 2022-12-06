function callback_GetSTRINGNetworkImg(src,~,genelist)

FigureHandle=src.Parent.Parent;
sce=guidata(FigureHandle);    

[speciestag] = i_selectspecies(3);