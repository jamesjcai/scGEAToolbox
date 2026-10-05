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

% Ask before any of the pickers open, the way the other annotation handlers
% do. It comes after the workspace check above, so a run that can only end in
% 'No SCE variables in Workspace.' does not ask about an overwrite that was
% never on offer. The stash below is what makes the import recoverable; this
% is about the active labels changing.
if ~gui.i_confirmoverwritecelltype(FigureHandle, sce), return; end

if ~allcell
    % What a separate 'Select a cell subtype, then an SCE variable that
    % contains the subtype annotation. Continue?' question used to say, now
    % said by the picker it was introducing. The picker's own Cancel is the
    % No that question was collecting, and the overwrite is already confirmed
    % above, so the extra round trip bought nothing.
    prompt = ['Select the cell subtype to re-annotate. The next step picks ' ...
        'the SCE that holds the annotation for those cells.'];

    celltypelist = natsort(unique(sce.c_cell_type_tx));
    if gui.i_isuifig(FigureHandle)
        % One type at a time: the line below compares C_CELL_TYPE_TX against
        % a single label, and the LISTDLG branch has always been 'single'.
        [indx, tf1] = gui.myListdlg(FigureHandle, celltypelist, ...
            'Select Cell Type', [], false, true, [], prompt);
    else
        [indx, tf1] = listdlg('PromptString', ...
            {'Select the cell subtype to re-annotate.', ...
             'The next step picks the SCE with its annotation.'}, ...
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
                'Select SCE:', [], false);
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
            'Select SCE Data File', 'MultiSelect', 'off');
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


% Said before the try. A data file is loaded without the size filter the
% workspace list applies, and a wrong-size one used to hit a bare ASSERT.
if insce.NumCells ~= sum(selecteidx)
    gui.myErrordlg(FigureHandle, sprintf(['The source has %d cells but ' ...
        '%d cells are to be labelled. Import needs an SCE holding the ' ...
        'same cells as this one.'], insce.NumCells, sum(selecteidx)));
    return;
end

try

    % A matching cell count is not a matching set of cells. When both sides
    % carry cell IDs, line the incoming labels up on those rather than on
    % column order, which is all the count above agrees on.
    [newtx, ok, idnote] = in_alignbyid(FigureHandle, sce, selecteidx, insce);
    if ~ok, return; end

    % Keep the labels this import replaces, the way the other annotation
    % handlers do (GUI.CALLBACK_ASSIGNCELLTYPEFROMATTRIB,
    % GUI.CALLBACK_SUBTYPEANNOTATION). The incoming SCE carries no trace of
    % what was here before, so without the stash an 'All Cells' import simply
    % erases the existing annotation. Taken after the size check, so an import
    % that is rejected leaves no stray 'old_cell_type_N' attribute behind.
    stashname = pkg.i_stashcelltypehistory(sce);

    sce.c_cell_type_tx(selecteidx) = newtx;
    gui.myGuidata(FigureHandle, sce, src);
    requirerefresh = true;
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

s = string(sprintf('Cell type annotation imported for %s.', ...
    pkg.i_plural(sum(selecteidx), 'cell')));
if strlength(idnote) > 0
    s = s + " " + idnote;
end
gui.myHelpdlg(FigureHandle, s + gui.i_stashnotice(stashname));

end % end of function

function [newtx, ok, note] = in_alignbyid(FigureHandle, sce, selecteidx, insce)
% Resolve the incoming labels onto the selected cells. OK is false when the
% user, told the IDs disagree, chose not to go ahead on position alone. NOTE
% is what the IDs showed, "" when there was nothing to check against, for the
% caller to add to what it reports.

newtx = insce.c_cell_type_tx;
ok = true;
note = "";

targetid = strings(0, 1);
if ~isempty(sce.c_cell_id)
    targetid = string(sce.c_cell_id(:));
    targetid = targetid(selecteidx);
end

[order, status, msg] = pkg.i_matchcellsbyid(targetid, insce.c_cell_id);

switch status
    case "inorder"
        % Say so. The whole point of the check is that the user can tell a
        % confirmed match from an unchecked one, and silence reads as the
        % latter.
        note = msg;
    case "reordered"
        newtx = newtx(order);
        note = msg;
    case {"mismatch", "duplicates"}
        % Not a refusal. Cell IDs pick up batch suffixes when SCEs are
        % merged and differ between readers, so IDs that do not line up are
        % sometimes the same cells spelled differently, and this import
        % worked on position alone before it ever looked at them. Say what
        % was found and default to No.
        answer = gui.myQuestdlg(FigureHandle, ...
            msg + " Import by position anyway?", ...
            'Cell IDs Do Not Match', {'Yes', 'No'}, 'No');
        ok = strcmp(answer, 'Yes');
    otherwise
        % "positional": one side has no IDs to check against, and this
        % import matched on column order alone before the check existed.
end
end
