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

% if ~gui.i_setwrkdir(preftagname), return; end
% s = getpref('scgeatoolbox', preftagname);
% s1 = sprintf('%s_workingfolder', extprogname);
% wkdir = fullfile(s, s1);
% 
% if ~exist(wkdir,"dir")
%     mkdir(wkdir);
% else
%     answer = gui.myQuestdlg(FigureHandle, 'Directory existing. Overwrite?');
%     if ~strcmp(answer,'Yes'), return; end
% end
% 
% fprintf('CURRENTWDIR = "%s"\n', wkdir);



a = findall(FigureHandle, 'type', 'axes');
h = findall(a, 'type', 'scatter');
ptsSelected = logical(h.BrushData.');
if ~any(ptsSelected)    
    % gui.myHelpdlg(FigureHandle, 'Please use the brush in the axes toolbar to select root cell(s).', '');
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

% [ndim] = gui.i_choose2d3d;
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

%[ok] = gui.i_confirmscript('Run Pseudotime Analysis (Monocle3)?', ...
%    extprogname, 'r');
%if ~ok, return; end

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

[y, idx] = ismember({'monocle3_pseudotime'}, ...
    sce.list_cell_attributes(1:2:end));
if y
    sce.list_cell_attributes{idx*2} = t_mono3;
else
    sce.list_cell_attributes = [sce.list_cell_attributes, ...
        {'monocle3_pseudotime', t_mono3}];
end

[y, idx] = ismember({'monocle3_qvalue'}, ...
    sce.list_gene_attributes(1:2:end));
if y
    sce.list_gene_attributes{idx*2} = q_mono3;
else
    sce.list_gene_attributes = [sce.list_gene_attributes, ...
        {'monocle3_qvalue', q_mono3}];
end

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


    function in_checkselected(~, ~)
        ptsSelected = logical(h.BrushData.');
        uiresume(FigureHandle);
    end
end

