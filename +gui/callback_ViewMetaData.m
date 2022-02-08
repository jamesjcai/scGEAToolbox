function callback_ViewMetaData(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    inputdlg('Data Info:','Metadata Viewer',[10 50],{char(sce.metadata)});
end
