function [needupdatesce] = callback_RunMonocle3(src, ~)

needupdatesce = false;
extprogname = 'R_monocle3';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);

% if ~gui.i_setwrkdir(preftagname), return; end
% s = getpref('scgeatoolbox', preftagname);
% s1 = sprintf('%s_workingfolder', extprogname);
% wkdir = fullfile(s, s1);
% 
% if ~exist(wkdir,"dir")
%     mkdir(wkdir);
% else
%     answer = questdlg('Directory existing. Overwrite?');
%     if ~strcmp(answer,'Yes'), return; end
% end
% 
% fprintf('CURRENTWDIR = "%s"\n', wkdir);

FigureHandle = src.Parent.Parent;
a = findall(FigureHandle, 'type', 'axes');
h = findall(a, 'type', 'scatter');
ptsSelected = logical(h.BrushData.');
if ~any(ptsSelected)
    helpdlg('Please use the brush in the axes toolbar to select root cell(s).', '');
    return;
end
idx = find(ptsSelected);

[ndim] = gui.i_choose2d3d;
if isempty(ndim), return; end

answer = questdlg('Keep temporary working files?');
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


sce = guidata(FigureHandle);
fw = gui.gui_waitbar;
try
    [t_mono3, s_mono3, ~, q_mono3] = run.r_monocle3(sce.X, idx, ndim, wkdir, isdebug);
catch ME
    gui.gui_waitbar(fw, true);
    errordlg(ME.message);
    return;
end
gui.gui_waitbar(fw);
if isempty(t_mono3) || length(t_mono3) ~= sce.NumCells
    errordlg('MONOCLE3 running time error.');
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
else
    sce.struct_cell_embeddings.('monocle3d') = s_mono3;
end

guidata(FigureHandle, sce);
needupdatesce = true;


if ~(ismcc || isdeployed)
    labels = {'Save pseudotime T to variable named:'};
    vars = {'Tmonocleout'};
    values = {t_mono3};
    msgfig = export2wsdlg(labels, vars, values);
    uiwait(msgfig)
else
    gui.i_exporttable(table(t_mono3), true, 'Tmonocleout', ...
        'MonocleResTable');
end


    % answer = questdlg('View Monocle3 Minimum Spanning Tree (MST)?');
    % switch answer
    %     case 'Yes'
    %         f=figure;
    %         gui.i_gscatter3(s_mono3, grp2idx(t_mono3));
    %         hold on
    %         if size(m_mono3,1) == 3
    %             plot3(m_mono3(1,:), m_mono3(2,:), m_mono3(3,:),'k-');
    %         else
    %             plot3(m_mono3(1,:), m_mono3(2,:), 'k-');
    %         end
    %         colorbar
    %         defaultToolbar = findall(f, 'tag','FigureToolBar');
    %         gui.add_3dcamera(defaultToolbar, 'Scores');
    %         hc = colorbar;
    %         hc.Label.String = 'Pseudotime';
    %     otherwise
    %         return;
    % end

end
