function callback_ViewMetaData(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
a = inputdlg('Data Info:', 'Metadata Viewer', [15, 80], {char(sce.metadata)});
if ~isempty(a)
    newmetadata = strtrim(string(a{1}));
    if ~isequal(newmetadata, sce.metadata)
        answer = questdlg('Save metadata?', '');
        if strcmp(answer, 'Yes')
            sce.metadata = newmetadata;
            guidata(FigureHandle, sce);
        end
    end
end
end
