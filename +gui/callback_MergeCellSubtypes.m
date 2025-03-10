function [requirerefresh, s] = callback_MergeCellSubtypes(src, ~, ...
    sourcetag, allcell)
if nargin < 4, allcell = false; end
if nargin < 3, sourcetag = 1; end
requirerefresh = false;
s = "";

switch sourcetag
    case 1
        a = evalin('base', 'whos');
        b = struct2cell(a);
        valididx = ismember(b(4, :), 'SingleCellExperiment');
        if sum(valididx) < 1
            warndlg('No SCE variables in Workspace.', '');
            return;
        end
end

[FigureHandle, sce, isui] = gui.gui_getfigsce(src);


if ~allcell
    answer = questdlg('Select a cell subtype, then an SCE variable that contains the subtype annotation. Continue?');
    if ~strcmp(answer, 'Yes'), return; end

    celltypelist = natsort(unique(sce.c_cell_type_tx));
    [indx, tf1] = listdlg('PromptString', ...
        {'Select Tissue Type(s):'}, ...
        'SelectionMode', 'single', ...
        'ListString', celltypelist, ...
        'ListSize', [220, 300]);
    if tf1 ~= 1, return; end
    selectedtype = celltypelist(indx);
    selecteidx = sce.c_cell_type_tx == selectedtype;
else
    selecteidx = true(sce.NumCells, 1);
end

switch sourcetag
    case 1
        % a=evalin('base','whos');
        % b=struct2cell(a);
        % valididx=ismember(b(4,:),'SingleCellExperiment');
        % if sum(valididx)<1
        %     warndlg('No SCE variables in Workspace.','');
        %     return;
        % end
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
            warndlg('No valid SCE variables in Workspace.', '');
            return;
        end

        [indx, tf] = listdlg('PromptString', {'Select SCE:'}, ...
            'liststring', b(1, :), ...
            'SelectionMode', 'single', ...
            'ListSize', [220, 300]);
        if tf ~= 1, return; end
        try
            insce = evalin('base', a(indx).name);
        catch ME
            errordlg(ME.message);
            return;
        end
    case 2
        [fname, pathname] = uigetfile({'*.mat', 'SCE Data File (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Select SCE Data File', 'MultiSelect', 'off');
        if isequal(fname, 0), return; end
        try
            scefile = fullfile(pathname, fname);
            x = load(scefile, 'sce');
            insce = x.sce;
        catch ME
            errordlg(ME.message);
            return;
        end
end % end of sourcetag


try
    assert(insce.NumCells == sum(selecteidx));
    sce.c_cell_type_tx(selecteidx) = ...
        insce.c_cell_type_tx;
    guidata(FigureHandle, sce);
    requirerefresh = true;
catch ME
    errordlg(ME.message);
    return;
end

end % end of function
