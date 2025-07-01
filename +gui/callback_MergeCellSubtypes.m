function [requirerefresh, s] = callback_MergeCellSubtypes(src, ~, ...
    sourcetag, allcell)
if nargin < 4, allcell = false; end
if nargin < 3, sourcetag = 1; end
requirerefresh = false;
s = "";

[FigureHandle, sce] = gui.gui_getfigsce(src);

switch sourcetag
    case 1
        a = evalin('base', 'whos');
        b = struct2cell(a);
        valididx = ismember(b(4, :), 'SingleCellExperiment');
        if sum(valididx) < 1
            gui.myWarndlg(FigureHandle,'No SCE variables in Workspace.');
            return;
        end
end


if ~allcell
    answer = gui.myQuestdlg(FigureHandle, ['Select a cell subtype, then ' ...
        'an SCE variable that contains the subtype annotation. Continue?']);
    if ~strcmp(answer, 'Yes'), return; end

    celltypelist = natsort(unique(sce.c_cell_type_tx));
    if gui.i_isuifig(FigureHandle)
        [indx, tf1] = gui.myListdlg(FigureHandle, celltypelist, ...
            'Select Cell Type(s):');
    else                        
        [indx, tf1] = listdlg('PromptString', ...
            {'Select Cell Type(s):'}, ...
            'SelectionMode', 'single', ...
            'ListString', celltypelist, ...
            'ListSize', [220, 300]);
    end
    if tf1 ~= 1, return; end
    selectedtype = celltypelist(indx);
    selecteidx = sce.c_cell_type_tx == selectedtype;
else
    selecteidx = true(sce.NumCells, 1);
end

switch sourcetag
    case 1
        b = b(:, valididx);
        a = a(valididx);

        valididx = false(length(a), 1);
        for k = 1:length(a)
            insce = evalin('base', a(k).name);
            if sum(selecteidx) == insce.NumCells
                valididx(k) = true;
            end
        end

        b = b(:, valididx);
        a = a(valididx);

        if isempty(a)
            gui.myWarndlg(FigureHandle, ...
                'No valid SCE variables in Workspace.');
            return;
        end

        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, b(1,:), ...
                'Select SCE:');
        else                            
            [indx, tf] = listdlg('PromptString', {'Select SCE:'}, ...
                'liststring', b(1, :), ...
                'SelectionMode', 'single', ...
                'ListSize', [220, 300]);
        end

        if tf ~= 1, return; end
        try
            insce = evalin('base', a(indx).name);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
    case 2
        [fname, pathname] = uigetfile({'*.mat', 'SCE Data File (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Select SCE Data File', ...
            'MultiSelect', 'off');
        if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end
        if isequal(fname, 0), return; end
        try
            scefile = fullfile(pathname, fname);
            x = load(scefile, 'sce');
            insce = x.sce;
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
end % end of sourcetag


try
    assert(insce.NumCells == sum(selecteidx));
    sce.c_cell_type_tx(selecteidx) = ...
        insce.c_cell_type_tx;
    gui.myGuidata(FigureHandle, sce, src);
    requirerefresh = true;
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

end % end of function
