function [requirerefresh] = callback_RenameBatchID(src, ~)

requirerefresh = false;


[FigureHandle, sce] = gui.gui_getfigsce(src);


if isempty(sce.c_batch_id)
    sce.c_batch_id = string(ones(sce.NumCells, 1));
end

if ~isstring(sce.c_batch_id)
    sce.c_batch_id = string(sce.c_batch_id);
end

answer = gui.myQuestdlg(FigureHandle, ...
    'Rename one batch ID, or edit all names in a table?', ...
    '', {'Select One', 'Edit Table', 'Cancel'}, 'Select One');
switch answer
    case 'Select One'
    case 'Edit Table'
        % Old and new names side by side; giving two batches the same new
        % name merges them.
        newvals = gui.i_relabeldlg(FigureHandle, sce.c_batch_id, ...
            'Rename Batch IDs');
        if ~isempty(newvals)
            sce.c_batch_id = newvals;
            requirerefresh = true;
        end
        gui.myGuidata(FigureHandle, sce, src);
        return;
    otherwise
        gui.myGuidata(FigureHandle, sce, src);
        return;
end

[ci, cLi] = pkg.i_grp2idxsorted(sce.c_batch_id);

% %---------
% [cLisorted,idx]=natsort(cLi);
% cisorted=ci;
% for k=1:length(idx), cisorted(ci==idx(k))=k; end
% ci=cisorted;
% cLi=cLisorted;
% %----------

if gui.i_isuifig(FigureHandle)
    [indxx, tfx] = gui.myListdlg(FigureHandle, cLi, ...
        'Select batch ID', [], false);
else
    [indxx, tfx] = listdlg('PromptString', ...
        {'Select batch ID'}, ...
        'SelectionMode', 'single', ...
        'ListString', string(cLi), 'ListSize', [220, 300]);
end
if tfx == 1
    i = ismember(ci, indxx);

    if gui.i_isuifig(FigureHandle) % Assuming FigureHandle is available
        newctype = gui.myInputdlg({'New batch ID name'}, 'Rename', ...
            cLi(ci(i)), FigureHandle);
    else
        newctype = inputdlg('New batch ID name', 'Rename', [1, 50], cLi(ci(i)));
    end


    if ~isempty(newctype)
        cLi(ci(i)) = newctype;
        sce.c_batch_id = string(cLi(ci));
        requirerefresh = true;
    end
end
gui.myGuidata(FigureHandle, sce, src);
end
