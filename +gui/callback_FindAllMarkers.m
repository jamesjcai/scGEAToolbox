function callback_FindAllMarkers(src, ~)

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

answer = gui.myQuestdlg(FigureHandle, 'Select Method', ...
'', {'Marker Gene Heatmap', 'Find All Markers'}, ...
'Marker Gene Heatmap');
switch answer
    case 'Find All Markers'
        in_findAllMarkers(sce, FigureHandle);
    case 'Marker Gene Heatmap'
        in_MarkerGeneHeatmap(sce, FigureHandle);
end
end


function in_findAllMarkers(sce, FigureHandle)

[thisc, ~] = gui.i_select1class(sce,[],[],[],FigureHandle);
if isempty(thisc), return; end
if isscalar(unique(thisc))
    gui.myWarndlg(FigureHandle, "All cells are in the same group.");
    return;
end
[T] = pkg.e_findallmarkers(sce.X, sce.g, thisc, [], [], [], true);
if ~isempty(T)
    % needwait = true;
    % [answer, filename] = gui.i_exporttable(T, needwait, 'Tallmarkers', ...
    %     'AllMarkersTable',[],[], FigureHandle);
    %         % "Tcellattrib","CellAttribTable"
    %         % "Tviolindata","ViolinPlotTable"
    %         % "Tcrosstabul","CrosstabulTable"
    %         % "Tcellsignmt","CellSignatTable"
    %         % "Tdpgenesres","DPGenesResTable"
    %         % "Tallmarkers","AllMarkersTable"
    % if ~isempty(answer)
    %     disp(filename);
    %     % gui.myHelpdlg(FigureHandle, sprintf('All Markers Table saved.'), '');
    %     if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'View marker table?'))
    %         gui.TableViewerApp(T, FigureHandle);
    %     end
    % end
    fw = gui.myWaitbar(FigureHandle);
    gui.TableViewerApp(T, FigureHandle);
    gui.myWaitbar(FigureHandle, fw);
else
    gui.myHelpdlg(FigureHandle, 'No results.', '');
end


end

function in_MarkerGeneHeatmap(sce, FigureHandle)
mfolder = fileparts(mfilename('fullpath'));
% unique(sce.c_cluster_id)

answer = gui.myQuestdlg(FigureHandle, "Only consider known (PangloaDB) marker genes?","");
if strcmp(answer, 'Yes')
    markergenes = pkg.i_get_panglaodbmarkers;
    idx = ismember(upper(sce.g), upper(markergenes));
    if sum(idx) > 2000
        fprintf('Size of input matrix: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);
        sce = sce.selectgenesbyindex(idx);
        sce = sce.qcfilter; % OK
        fprintf('Size of filtered matrix: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);
    else
        if ~strcmp('Yes', gui.myQuestdlg(FigureHandle,...
            "Not enough marker genes (n < 2000). Use all genes to search for markers."))
            return;
        end
    end
elseif strcmp(answer, 'No')
    disp('Consider all genes for marker gene search.');
else
    return;
end

[thisc, clabel] = gui.i_select1class(sce,[],[],[],FigureHandle);
if isempty(thisc), return; end
if isscalar(unique(thisc))
    gui.myWarndlg(FigureHandle, "All cells are in the same group.");
    return;
end
% [c, cL, noanswer] = gui.i_reordergroups(thisc, [], FigureHandle);
% if noanswer, return; end

% Which cell types go on the heatmap. Every level of the grouping variable
% used to go on it, and an annotation with thirty types -- or with a
% handful of cells sitting in a type nobody asked about -- gave a heatmap
% too dense to read, with ten marker genes picked per type either way.
% The order the list is offered in is the user's: a long annotation is
% easiest to work through by size, a familiar one alphabetically.
[picked, levels] = gui.i_selectgroupsubset(thisc, clabel, FigureHandle, "ask");
if isempty(picked), return; end
if numel(levels) < 2
    gui.myWarndlg(FigureHandle, "Select at least two cell types.");
    return;
end
if ~all(picked)
    numlevels = numel(unique(string(thisc)));
    sce = sce.selectcells(picked);
    thisc = thisc(picked);
    fprintf('%d of %d groups kept (%d of %d cells).\n', ...
        numel(levels), numlevels, sce.NumCells, numel(picked));
end

[c] = findgroups(thisc);
answer = gui.myQuestdlg(FigureHandle, 'Generate marker gene heatmap', ...
'Select Method', {'Method 1 (DE 🐇)', 'Method 2 (scGeneFit 🐢)', ...
'Method 3 (LASSO 🐢🐢)'}, 'Method 1 (DE 🐇)');
switch answer
    case 'Method 1 (DE 🐇)'
        methodid = 1;
    case 'Method 2 (scGeneFit 🐢)'
        methodid = 2;
    case 'Method 3 (LASSO 🐢🐢)'
        methodid = 3;
    otherwise
        return;
end

fw = gui.myWaitbar(FigureHandle);

% speciestag = gui.i_selectspecies(2, false, FigureHandle);
% if isempty(speciestag)
%     requirerefresh = false;
%     return;
% end
speciestag = "human";
load(fullfile(mfolder, ...
'..', 'assets', 'Biomart', sprintf('Biomart_%s_genes.mat',speciestag)), 'T');
ApprovedSymbol = string(T.GeneName);
sce = sce.rmmtgenes;
sce = sce.rmhemoglobingenes;
sce = sce.rmribosomalgenes;
sce = sce.rmlncrnagenes;
[idx] = ~ismember(upper(sce.g), upper(ApprovedSymbol));
if any(idx)
    sce = sce.selectgenesbyindex(~idx);
end

fprintf('Size of matrix used for search: %d genes x %d cells\n', ...
sce.NumGenes, sce.NumCells);

try
    [markerlist] = sc_pickmarkers(sce.X, sce.g, c, 10, methodid);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

glist = string(markerlist{1}(:));
for k = 2:length(markerlist)
    glist = [glist; string(markerlist{k}(:))];
end
gui.myWaitbar(FigureHandle, fw);
gui.i_heatmap(sce, glist, thisc, FigureHandle);
end
