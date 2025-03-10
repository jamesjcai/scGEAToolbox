function [requirerefresh] = callback_RenameBatchID(src, ~)

requirerefresh = false;

[FigureHandle, sce, isui] = gui.gui_getfigsce(src);

if isempty(sce.c_batch_id)
    sce.c_batch_id = string(ones(sce.NumCells, 1));
    %errordlg('sce.c_batch_id undefined');
    %return;
end
%answer = gui.myQuestdlg(FigureHandle, 'Rename batch ID?');
%if ~strcmp(answer, 'Yes'), return; end

if ~isstring(sce.c_batch_id)
    sce.c_batch_id = string(sce.c_batch_id);
end
[ci, cLi] = pkg.i_grp2idxsorted(sce.c_batch_id);

% %---------
% [cLisorted,idx]=natsort(cLi);
% cisorted=ci;
% for k=1:length(idx), cisorted(ci==idx(k))=k; end
% ci=cisorted;
% cLi=cLisorted;
% %----------


[indxx, tfx] = listdlg('PromptString', ...
    {'Select batch ID'}, ...
    'SelectionMode', 'single', ...
    'ListString', string(cLi), 'ListSize', [220, 300]);
if tfx == 1
    i = ismember(ci, indxx);
    newctype = inputdlg('New batch ID name', 'Rename', ...
        [1, 50], cLi(ci(i)));
    if ~isempty(newctype)
        cLi(ci(i)) = newctype;
        sce.c_batch_id = string(cLi(ci));
        requirerefresh = true;
    end
end
guidata(FigureHandle, sce);
end
