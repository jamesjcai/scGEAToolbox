function callback_ViewMetaData(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);

if gui.i_isuifig(FigureHandle)
    % a = gui.myInputdlg({'Data Info:'}, 'Metadata Viewer', {char(sce.metadata)}, FigureHandle);
    a = gui.myTextareadlg(FigureHandle, {'Data Info:'}, 'Metadata Viewer', {sce.metadata}, true);
else
    a = inputdlg('Data Info:', 'Metadata Viewer', [15, 80], {char(sce.metadata)});
end


if ~isempty(a)
    newmetadata = strtrim(string(a{1}));
    if ~isequal(newmetadata, sce.metadata)
        answer = gui.myQuestdlg(FigureHandle, 'Save metadata?', '');
        if strcmp(answer, 'Yes')
            sce.metadata = newmetadata;
            guidata(FigureHandle, sce);
        end
    end
end
end
