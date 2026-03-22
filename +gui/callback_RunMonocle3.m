function [needupdatesce] = callback_RunMonocle3(src, ~)

needupdatesce = false;
extprogname = 'R_monocle3';
[FigureHandle, sce] = gui.gui_getfigsce(src);

[wrkdir] = gui.i_getwrkdir;
if isempty(wrkdir)
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end
end


a = findall(FigureHandle, 'type', 'axes');
h = findall(a, 'type', 'scatter');
ptsSelected = logical(h.BrushData.');
if ~any(ptsSelected)
    answer = gui.myQuestdlg(FigureHandle, 'Use brush to select root cell(s). Ready?','');
    if ~strcmp(answer, 'Yes'), return; end
    b = brush(FigureHandle);
    b.ActionPostCallback=@in_checkselected;
    b.Enable="on";
    disp('Waiting for user to finish brushing...');
    uiwait(FigureHandle);  % Pauses the execution until uiresume is called
    disp('User finished brushing!');
    b.Enable="off";
    if any(ptsSelected)
        answer = gui.myQuestdlg(FigureHandle, 'Root cell(s) selected. Continue?','');
        if ~strcmp(answer, 'Yes'), return; end
    else
        gui.myHelpdlg(FigureHandle, 'No root cell(s) are selected.');
        return;
    end
end
idx = find(ptsSelected);
if isempty(idx), gui.myWarndlg(FigureHandle, 'Root cell(s) is missing.'); return; end

ndim = 2;
if isempty(ndim), return; end

answer = gui.myQuestdlg(FigureHandle, 'Keep temporary working files?');
switch answer
    case 'Yes'
        isdebug = true;
    case 'No'
        isdebug = false;
    case 'Cancel'
        return;
    otherwise
        return;
end

fw = gui.myWaitbar(FigureHandle);
try
    [t_mono3, s_mono3, ~, q_mono3] = run.r_monocle3(sce.X, idx, ndim, wrkdir, isdebug);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
if isempty(t_mono3) || length(t_mono3) ~= sce.NumCells
    gui.myErrordlg(FigureHandle, 'MONOCLE3 running time error.');
    return;
end

sce.setCellAttribute('monocle3_pseudotime', t_mono3);

sce.setCellAttribute('monocle3_qvalue', q_mono3);

if ndim == 2
    sce.struct_cell_embeddings.('monocle2d') = s_mono3;
elseif ndim == 3
    sce.struct_cell_embeddings.('monocle3d') = s_mono3;
else
    error('Invalid ndim. (Valid options: 2 or 3)');
end

gui.myGuidata(FigureHandle, sce, src);
needupdatesce = true;

gui.myHelpdlg(FigureHandle, 'Monocle3 pseudotime T and embedding S have been saved in SCE.');

[y, idx] = ismember({'monocle3_pseudotime'}, sce.list_cell_attributes(1:2:end));
if y
    answer = gui.myQuestdlg(FigureHandle, ...
        ['Color cells using pseudotime T and show' ...
        ' Monocle3 embedding S?'], ...
        '',{'Yes','Color cells only','Show embedding only'},'Yes');
    switch answer
        case 'Yes'
            sce.c = sce.list_cell_attributes{idx*2};
            sce.s = sce.struct_cell_embeddings.('monocle2d');
        case 'Color cells only'
            sce.c = sce.list_cell_attributes{idx*2};
        case 'Show embedding only'
            sce.s = sce.struct_cell_embeddings.('monocle2d');
        otherwise
    end
    gui.myGuidata(FigureHandle, sce, src);
end


function in_checkselected(~, ~)
        ptsSelected = logical(h.BrushData.');
        uiresume(FigureHandle);
    end
end
