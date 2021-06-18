function varargout = sc_scatter_sce(sce, varargin)

if nargin < 1
    error('Usage: sc_scatter_sce(sce)');
end
if ~isa(sce, 'SingleCellExperiment')
    error('requires sce=SingleCellExperiment();');
end

import pkg.*
import gui.*

p = inputParser;
checkCS = @(x) isempty(x) | size(sce.X, 2) == length(x);
addRequired(p, 'sce', @(x) isa(x, 'SingleCellExperiment'));
addOptional(p, 'c', sce.c, checkCS);
addOptional(p, 's', [], checkCS);
addOptional(p, 'methodid', 1, @isnumeric);
parse(p, sce, varargin{:});
cin = p.Results.c;
sin = p.Results.s;
methodid = p.Results.methodid;
ax = [];
bx = [];

if isempty(cin)
    sce.c = ones(size(sce.X, 2), 1);
else
    sce.c = cin;
end
if ~isempty(sin)
    sce.s = sin;
end

[c, cL] = grp2idx(sce.c);

FigureHandle = figure('Name', 'SC_SCATTER', ...
    'position', round(1.5 * [0 0 560 420]), ...
    'visible', 'off');
movegui(FigureHandle, 'center');

hAx = axes('Parent', FigureHandle);
[h] = gui.i_gscatter3(sce.s, c, methodid);
title(sce.title);

dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

% defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle
defaultToolbar = findall(FigureHandle, 'Type', 'uitoolbar');

% UitoolbarHandle2 = uitoolbar( 'Parent', FigureHandle ) ;
% set( UitoolbarHandle2, 'Tag' , 'FigureToolBar2' , ...
%     'HandleVisibility' , 'on' , ...
%     'Visible' , 'on' ) ;

UitoolbarHandle = uitoolbar('Parent', FigureHandle);
set(UitoolbarHandle, 'Tag', 'FigureToolBar', ...
    'HandleVisibility', 'off', ...
    'Visible', 'on');

mfolder = fileparts(mfilename('fullpath'));

% UitoolbarHandle = uitoolbar(FigureHandle);
pt3 = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'list.gif'));
ptImage = ind2rgb(img, map);
pt3.CData = ptImage;
pt3.Tooltip = 'Select a gene to show expression';
pt3.ClickedCallback = @callback_ShowGeneExpr;

pt3a = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'list2.gif'));
ptImage = ind2rgb(img, map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Show cell states';
pt3a.ClickedCallback = @ShowCellStats;

pt3a = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-pointfig.gif'));
ptImage = ind2rgb(img, map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Select cells by class';
pt3a.ClickedCallback = @callback_SelectCellsByClass;

pt3a = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-effects.gif'));
ptImage = ind2rgb(img, map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Filter genes and cells';
pt3a.ClickedCallback = @SelectCellsByQC;

% ------------------

ptlabelclusters = uitoggletool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(matlabroot, ...
    'toolbox', 'matlab', 'icons', 'plotpicker-scatter.gif'));
% map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
ptImage = ind2rgb(img, map);
ptlabelclusters.CData = ptImage;
ptlabelclusters.Tooltip = 'Label clusters';
ptlabelclusters.ClickedCallback = @LabelClusters;

% ------------------ clustering

ptaddcluster = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-glyplot-face.gif'));
ptImage = ind2rgb(img, map);
ptaddcluster.CData = ptImage;
ptaddcluster.Tooltip = 'Add brushed cells to a new cluster';
ptaddcluster.ClickedCallback = @Brushed2NewCluster;

ptmergecluster = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-pzmap.gif'));
ptImage = ind2rgb(img, map);
ptmergecluster.CData = ptImage;
ptmergecluster.Tooltip = 'Merge brushed cells to same cluster';
ptmergecluster.ClickedCallback = @Brushed2MergeClusters;

ptShowClu = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-geoscatter.gif'));
ptImage = ind2rgb(img, map);
ptShowClu.CData = ptImage;
ptShowClu.Tooltip = 'Show clusters individually';
ptShowClu.ClickedCallback = @ShowClustersPop;

ptcluster = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-dendrogram.gif'));
ptImage = ind2rgb(img, map);
ptcluster.CData = ptImage;
ptcluster.Tooltip = 'Clustering using embedding S';
ptcluster.ClickedCallback = @ClusterCellsS;

ptcluster = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-gscatter.gif'));
ptImage = ind2rgb(img, map);
ptcluster.CData = ptImage;
ptcluster.Tooltip = 'Clustering using expression matrix X';
ptcluster.ClickedCallback = @ClusterCellsX;

% -------------

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'brush.gif'));
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Cell types of brushed cells';
pt5.ClickedCallback = @Brush4Celltypes;

ptclustertype = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(matlabroot, ...
    'toolbox', 'matlab', 'icons', 'plotpicker-contour.gif'));
ptImage = ind2rgb(img, map);
ptclustertype.CData = ptImage;
ptclustertype.Tooltip = 'Cell types of clusters';
ptclustertype.ClickedCallback = @DetermineCellTypeClusters;

pt4 = uipushtool(UitoolbarHandle, 'Separator', 'off');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-scatterhist.gif'));
ptImage = ind2rgb(img, map);
pt4.CData = ptImage;
pt4.Tooltip = 'Rename cell type';
pt4.ClickedCallback = @RenameCellType;

pt4 = uipushtool(UitoolbarHandle, 'Separator', 'off');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-kagi.gif'));
ptImage = ind2rgb(img, map);
pt4.CData = ptImage;
pt4.Tooltip = 'Marker genes of brushed cells';
pt4.ClickedCallback = @callback_Brush4Markers;

pt4mrkheat = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-plotmatrix.gif'));
ptImage = ind2rgb(img, map);
pt4mrkheat.CData = ptImage;
pt4mrkheat.Tooltip = 'Marker gene heatmap';
pt4mrkheat.ClickedCallback = @callback_MarkerGeneHeatmap;

% --------------------------

ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-candle.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Compare Differentiation Potency';
ptpseudotime.ClickedCallback = @callback_ComparePotency;

ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-arxtimeseries.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Run pseudotime analysis (Monocle)';
ptpseudotime.ClickedCallback = @RunTrajectoryAnalysis;

ptpseudotime = uipushtool(defaultToolbar, ...
    'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-comet.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Plot pseudotime trajectory';
ptpseudotime.ClickedCallback = @DrawTrajectory;

ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-priceandvol.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Compare Gene Expression between Classes';
ptpseudotime.ClickedCallback = @callback_CompareGeneBtwCls;

pt4 = uipushtool(defaultToolbar, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-boxplot.gif'));
ptImage = ind2rgb(img, map);
pt4.CData = ptImage;
pt4.Tooltip = 'Compare 2 groups (DE analysis)';
pt4.ClickedCallback = @callback_DEGene2Groups;

ptpseudotime = uipushtool(defaultToolbar, ...
    'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-andrewsplot.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Function enrichment of HVG genes';
ptpseudotime.ClickedCallback = @callback_GSEA_HVGs;

ptnetwork = uipushtool(defaultToolbar, ...
    'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'noun_Network_691907.gif'));
ptImage = ind2rgb(img, map);
ptnetwork.CData = ptImage;
ptnetwork.Tooltip = 'Build gene regulatory network';
ptnetwork.ClickedCallback = @callback_BuildGeneNetwork;

ptnetwork = uipushtool(defaultToolbar, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'noun_Deep_Learning_2424485.gif'));
ptImage = ind2rgb(img, map);
ptnetwork.CData = ptImage;
ptnetwork.Tooltip = 'Compare two scGRNs';
ptnetwork.ClickedCallback = @callback_CompareGeneNetwork;

pt2 = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-qqplot.gif'));
ptImage = ind2rgb(img, map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delete selected cells';
pt2.ClickedCallback = @DeleteSelectedCells;

pt = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, 'resources', 'export.gif'));
ptImage = ind2rgb(img, map);
pt.CData = ptImage;
pt.Tooltip = 'Export & save data';
pt.ClickedCallback = @callback_SaveX;

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-geobubble.gif'));
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Embedding';
pt5.ClickedCallback = @EmbeddingAgain;

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'multiscale.gif'));
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Multi-embedding view';
pt5.ClickedCallback = @gui.callback_MultiEmbeddings;

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-image.gif'));  % plotpicker-pie
% map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Switch 2D/3D';
pt5.ClickedCallback = @Switch2D3D;

pt5pickmk = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-rose.gif'));  % plotpicker-pie
ptImage = ind2rgb(img, map);
pt5pickmk.CData = ptImage;
pt5pickmk.Tooltip = 'Switch scatter plot marker type';
pt5pickmk.ClickedCallback = @callback_PickPlotMarker;

pt5pickcl = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-compass.gif'));  % plotpicker-pie
ptImage = ind2rgb(img, map);
pt5pickcl.CData = ptImage;
pt5pickcl.Tooltip = 'Switch color maps';
pt5pickcl.ClickedCallback = {@gui.callback_PickColorMap, ...
    numel(unique(c))};

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'plotpicker-geobubble2.gif'));
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Refresh';
pt5.ClickedCallback = @RefreshAll;

gui.add_3dcamera(defaultToolbar, 'AllCells');



m = uimenu(FigureHandle,'Text','Experimental');
uimenu(m,'Text','Remove batch effect using Harmony (python required)...','Callback',@HarmonyPy);
uimenu(m,'Text','Extract cells by marker(+/-) expression...',...
    'Callback',@callback_SelectCellsByMarker);
uimenu(m,'Text','Ligand-receptor mediated intercellular crosstalk...',...
    'Callback',@callback_DetectIntercellularCrosstalk);
uimenu(m,'Text','Doublet Detection (python required)...',...
    'Callback',@DoubletDetection);
uimenu(m,'Text','Gene Expression Statistics...',...
    'Callback',@callback_CalculateGeneStats);

% handles = guihandles( FigureHandle ) ;
% guidata( FigureHandle, handles ) ;
set(FigureHandle, 'visible', 'on');
guidata(FigureHandle, sce);

% set(FigureHandle,'CloseRequestFcn',@closeRequest);

if nargout > 0
    varargout{1} = FigureHandle;
end

% ------------------------
% Callback Functions
% ------------------------

% function closeRequest(hObject,~)
% ButtonName = questdlg('Close SC_SCATTER?', ...
%                          '', ...
%                          'Yes','No','No');
% switch ButtonName
%     case 'Yes'
%         delete(hObject);
%     case 'No'
%         return;
% end
% end

    function SelectCellsByQC(src, ~)
        [requirerefresh,highlightindex]=gui.callback_SelectCellsByQC(src);
        if requirerefresh
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c);
            RefreshAll(src, 1, true);            
        end
        if ~isempty(highlightindex)
            h.BrushData=highlightindex;
            % set(h,'BrushData',highlightindex');
        end
    end

    function HarmonyPy(src, ~)
        if gui.callback_Harmonypy(src)
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c);
            RefreshAll(src, 1, true, false);
            ButtonName = questdlg('Update Saved Embedding?', ...
                '', ...
                'tSNE','UMAP','PHATE','tSNE');
            methodtag=lower(ButtonName);
            if ismember(methodtag,{'tsne','umap','phate'})
                sce.struct_cell_embeddings.(methodtag)=sce.s;
            end
        end
    end


    function DoubletDetection(src, ~)
        
        [isDoublet,doubletscore,done]=gui.callback_DoubletDetection(src);
        sce.list_cell_attributes=[sce.list_cell_attributes,...
                {'doublet_score',doubletscore}];
        guidata(FigureHandle,sce);
        if done && any(isDoublet)
            answer=questdlg(sprintf('Delete detected doublets (n=%d)?',...
                sum(isDoublet)),...
                '','Yes','No, show scores','Cancel','Yes');
            switch answer
                case 'Yes'
                    i_deletecells(isDoublet);
                case 'No'
                    i_showstate(doubletscore);
                    return;
                case 'Cancel'
                    return;
                otherwise
                    retrun;
            end
            %sce = guidata(FigureHandle);
            %sce.c=isDoublet;
            %[c, cL] = grp2idx(sce.c);
            %RefreshAll(src, 1, true, false);
        elseif done && ~any(isDoublet)
            helpdlg('No doublet found.');
            %sce.c=doubletscore;
            %RefreshAll(src, 1, true, false);
        end

    end

% =========================
    function RefreshAll(src, ~, keepview, keepcolr)
        if nargin < 4
            keepcolr = false;
        end
        if nargin < 3
            keepview = false;
        end
        
        if keepview || keepcolr
            [para] = i_getoldsettings(src);
        end
        if size(sce.s, 2) > 2
            if ~isempty(h.ZData)
                if keepview
                    [ax, bx] = view();
                end
                h = gui.i_gscatter3(sce.s, c, methodid);
                if keepview
                    view(ax, bx);
                end
            else
                h = gui.i_gscatter3(sce.s(:, 1:2), c, methodid);
            end
        else
            h = gui.i_gscatter3(sce.s(:, 1:2), c, methodid);
        end
        if keepview
            h.Marker = para.oldMarker;
            h.SizeData = para.oldSizeData;
        end
        if keepcolr
            colormap(para.oldColorMap);
        else
            kc = numel(unique(c));
            if kc <= 50
                colormap(lines(kc));
            else
                colormap default;
            end
        end
        title(sce.title);
        pt5pickcl.ClickedCallback = {@callback_PickColorMap, ...
            numel(unique(c))};
        guidata(FigureHandle, sce);
        ptlabelclusters.State = 'off';
        % UitoolbarHandle.Visible='off';
        % UitoolbarHandle.Visible='on';
    end

    function Switch2D3D(src, ~)
        [para] = i_getoldsettings(src);
        if isempty(h.ZData)   % current 2 D
            if ~(size(sce.s, 2) > 2)
                helpdlg('Canno swith to 3-D. SCE.S is 2-D');
                return
            end
            h = gui.i_gscatter3(sce.s, c, methodid);
            if ~isempty(ax) && ~isempty(bx) && ~any([ax bx] == 0)
                view(ax, bx);
            else
                view(3);
            end
        else                 % current 3D do following
            [ax, bx] = view();
            answer = questdlg('Which view to be used to project cells?', '', ...
                'Default View', 'Current View', 'Cancel', 'Default View');
            switch answer
                case 'Default View'
                    h = gui.i_gscatter3(sce.s(:, 1:2), c, methodid);
                case 'Current View'
                    sx = pkg.i_3d2d(sce.s, ax, bx);
                    h = gui.i_gscatter3(sx(:, 1:2), c, methodid);
                case {'Cancel', ''}
                    return
            end
        end
        title(sce.title);
        h.Marker = para.oldMarker;
        h.SizeData = para.oldSizeData;
        colormap(para.oldColorMap);
    end

    function RenameCellType(~, ~)
        if isempty(sce.c_cell_type_tx)
            errordlg('sce.c_cell_type_tx undefined');
            return
        end
        answer = questdlg('Rename a cell type?');
        if ~strcmp(answer, 'Yes')
            return
        end
        [ci, cLi] = grp2idx(sce.c_cell_type_tx);
        [indxx, tfx] = listdlg('PromptString', {'Select a cell type', ...
            '', ''}, 'SelectionMode', 'single', 'ListString', string(cLi));
        if tfx == 1
            i = ismember(ci, indxx);
            newctype = inputdlg('New cell type', 'Rename', [1 50], cLi(ci(i)));
            if ~isempty(newctype)
                cLi(ci(i)) = newctype;
                sce.c_cell_type_tx = cLi(ci);
                [c, cL] = grp2idx(sce.c_cell_type_tx);
                i_labelclusters(false);
            end
        end
        guidata(FigureHandle, sce);
    end

    function EmbeddingAgain(src, ~)
        answer = questdlg('Which embedding method?', 'Select method', 'tSNE', 'UMAP', 'PHATE', 'tSNE');
        if ~ismember(answer, {'tSNE', 'UMAP', 'PHATE'})
            return
        end
        if isempty(sce.struct_cell_embeddings)
            sce.struct_cell_embeddings = struct('tsne', [], 'umap', [], 'phate', []);
        end
        methodtag = lower(answer);
        usingold = false;
        if ~isempty(sce.struct_cell_embeddings.(methodtag))
            answer1 = questdlg(sprintf('Use existing %s embedding or re-compute new embedding?', ...
                upper(methodtag)), '', ...
                'Use existing', 'Re-compute', 'Cancel', 'Use existing');
            switch answer1
                case 'Use existing'
                    sce.s = sce.struct_cell_embeddings.(methodtag);
                    usingold = true;
                case 'Re-compute'
                    usingold = false;
                case {'Cancel', ''}
                    return
            end
        end
        if ~usingold
            answer2 = questdlg(sprintf('Use highly variable genes (HVGs, n=2000) or use all genes (n=%d)?', sce.NumGenes), ...
                '', ...
                'Top 2000 HVGs', 'All Genes', 'Cancel', 'Top 2000 HVGs');
            switch answer2
                case 'All Genes'
                    usehvgs = false;
                case 'Top 2000 HVGs'
                    usehvgs = true;
                case {'Cancel', ''}
                    return
            end
            fw = gui.gui_waitbar;
            try
                forced = true;
                sce = sce.embedcells(methodtag, forced, usehvgs);
            catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
                return
            end
            gui.gui_waitbar(fw);
        end
        RefreshAll(src, 1, true, false);
        guidata(FigureHandle, sce);
    end

    function DetermineCellTypeClusters(~, ~)
        answer = questdlg('Assign cell type identity to clusters?');
        if ~strcmp(answer, 'Yes')
            return
        end
        
        answer = questdlg('Which species?', 'Select Species', 'Mouse', 'Human', 'Mouse');
        
        if strcmp(answer, 'Human')
            speciestag = "human";
        elseif strcmp(answer, 'Mouse')
            speciestag = "mouse";
        else
            return
        end
        organtag = "all";
        databasetag = "panglaodb";
        dtp = findobj(h, 'Type', 'datatip');
        delete(dtp);
        cLdisp = cL;
        for i = 1:max(c)
            ptsSelected = c == i;
            [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, ...
                sce.s, ptsSelected, ...
                speciestag, organtag, databasetag);
            ctxt = Tct.C1_Cell_Type;
            
            [indx, tf] = listdlg('PromptString', {'Select cell type', ...
                '', ''}, 'SelectionMode', 'single', 'ListString', ctxt);
            if tf == 1
                ctxt = Tct.C1_Cell_Type{indx};
            else
                return
            end
            hold on;
            ctxtdisp = strrep(ctxt, '_', '\_');
            ctxtdisp = sprintf('%s_{%d}', ctxtdisp, i);
            cLdisp{i} = ctxtdisp;
            
            ctxt = sprintf('%s_{%d}', ctxt, i);
            cL{i} = ctxt;
            
            row = dataTipTextRow('', cLdisp(c));
            h.DataTipTemplate.DataTipRows = row;
            if size(sce.s, 2) >= 2
                siv = sce.s(ptsSelected, :);
                si = mean(siv, 1);
                idx = find(ptsSelected);
                [k] = dsearchn(siv, si);
                datatip(h, 'DataIndex', idx(k));
                % text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
                %     'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
                %     elseif size(sce.s,2)==2
                %             si=mean(sce.s(ptsSelected,:));
                %             text(si(:,1),si(:,2),sprintf('%s',ctxt),...
                %                  'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            end
            hold off;
        end
        sce.c_cell_type_tx = string(cL(c));
        guidata(FigureHandle, sce);
    end

    function Brushed2NewCluster(~, ~)
        answer = questdlg('Make a new cluster out of brushed cells?');
        if ~strcmp(answer, 'Yes')
            return
        end
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return
        end
        c(ptsSelected) = max(c) + 1;
        [c, cL] = grp2idx(c);
        sce.c = c;
        [ax, bx] = view();
        [h] = gui.i_gscatter3(sce.s, c, methodid);
        title(sce.title);
        view(ax, bx);
        i_labelclusters(true);
        sce.c_cluster_id = c;
        guidata(FigureHandle, sce);
    end

    function Brushed2MergeClusters(~, ~)
        answer = questdlg('Merge brushed cells into one cluster?');
        if ~strcmp(answer, 'Yes')
            return
        end
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are brushed");
            return
        end
        c_members = unique(c(ptsSelected));
        if numel(c_members) == 1
            warndlg("All brushed cells are in one cluster");
            return
        else
            [indx, tf] = listdlg('PromptString', {'Select target cluster', ...
                '', ''}, 'SelectionMode', 'single', 'ListString', string(c_members));
            if tf == 1
                c_target = c_members(indx);
            else
                return
            end
        end
        c(ismember(c, c_members)) = c_target;
        [c, cL] = grp2idx(c);
        sce.c = c;
        [ax, bx] = view();
        [h] = gui.i_gscatter3(sce.s, c, methodid);
        title(sce.title);
        view(ax, bx);
        i_labelclusters(true);
        sce.c_cluster_id = c;
        guidata(FigureHandle, sce);
    end

    function Brush4Celltypes(~, ~)
        answer = questdlg('Label cell type of brushed cells?');
        if ~strcmp(answer, 'Yes')
            return
        end
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return
        end
        answer = questdlg('Which species?', 'Select Species', 'Mouse', 'Human', 'Mouse');
        switch answer
            case 'Human'
                speciestag = "human";
            case 'Mouse'
                speciestag = "mouse";
            otherwise
                return
        end
        fw = gui.gui_waitbar;
        [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, sce.s, ptsSelected, ...
            speciestag, "all", "panglaodb");
        ctxt = Tct.C1_Cell_Type;
        gui.gui_waitbar(fw);
        
        [indx, tf] = listdlg('PromptString', {'Select cell type', ...
            '', ''}, 'SelectionMode', 'single', 'ListString', ctxt);
        if tf == 1
            ctxt = Tct.C1_Cell_Type{indx};
        else
            return
        end
        
        hold on;
        ctxt = strrep(ctxt, '_', '\_');
        if size(sce.s, 2) >= 3
            scatter3(sce.s(ptsSelected, 1), sce.s(ptsSelected, 2), sce.s(ptsSelected, 3), 'x');
            si = mean(sce.s(ptsSelected, :));
            text(si(:, 1), si(:, 2), si(:, 3), sprintf('%s', ctxt), ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        elseif size(sce.s, 2) == 2
            scatter(sce.s(ptsSelected, 1), sce.s(ptsSelected, 2), 'x');
            si = mean(sce.s(ptsSelected, :));
            text(si(:, 1), si(:, 2), sprintf('%s', ctxt), ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        end
        hold off;
    end

    function ShowCellStats(~, ~)
        % FigureHandle=src.Parent.Parent;
        sce = guidata(FigureHandle);
        listitems = {'Library Size', 'Mt-reads Ratio', ...
            'Mt-genes Expression', 'HgB-genes Expression', ...
            'Cell Cycle Phase', ...
            'Cell Type', 'Cluster ID', 'Batch ID'};
        % if ~ismember('cell potency',sce.list_cell_attributes)
        %    listitems{end+1}='Cell Potency';
        % end
        for k = 1:2:length(sce.list_cell_attributes)
            listitems = [listitems, sce.list_cell_attributes{k}];
        end
        [indx, tf] = listdlg('PromptString', {'Select statistics', ...
            '', ''}, 'SelectionMode', 'single', 'ListString', listitems);
        if tf ~= 1
            return
        end
        switch indx
            case 1
                ci = sum(sce.X);
                ttxt = "Library Size";
                pkg.i_stem3scatter(sce.s(:, 1), sce.s(:, 2), ci, ttxt);
                view(2);
                return
            case 2
                fw = gui.gui_waitbar;
                i = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
                lbsz = sum(sce.X, 1);
                lbsz_mt = sum(sce.X(i, :), 1);
                ci = lbsz_mt ./ lbsz;
                ttxt = "mtDNA%";
                pkg.i_stem3scatter(sce.s(:, 1), sce.s(:, 2), ci, ttxt);
                view(2);
                colorbar;
                pause(1);
                gui.gui_waitbar(fw);
                return
            case 3
                idx = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
                n = sum(idx);
                if n > 0
                    [ax, bx] = view();
                    if n <= 9
                        i_markergenespanel(sce.X, sce.g, sce.s, ...
                            sce.g(idx), [], 9, ax, bx, 'Mt-genes');
                    else
                        i_markergenespanel(sce.X, sce.g, sce.s, ...
                            sce.g(idx), [], 16, ax, bx, 'Mt-genes');
                    end
                else
                    warndlg('No mt-genes found');
                end
                return
            case 4 % HgB-genes
                idx = startsWith(sce.g, 'hba-', 'IgnoreCase', true) | startsWith(sce.g, 'hbb-', 'IgnoreCase', true);
                if any(idx)
                    ttxt = sprintf("%s+", sce.g(idx));
                    ci = sum(sce.X(idx, :), 1);
                    pkg.i_stem3scatter(sce.s(:, 1), sce.s(:, 2), ci, ttxt);
                else
                    warndlg('No HgB-genes found');
                end
                return
            case 5   % "Cell Cycle Phase";
                if isempty(sce.c_cell_cycle_tx)
                    fw = gui.gui_waitbar;
                    sce = sce.estimatecellcycle;
                    gui.gui_waitbar(fw);
                end
                [ci, tx] = grp2idx(sce.c_cell_cycle_tx);
                ttxt = sprintf('%s|', string(tx));
            case 6 % cell type
                ci = sce.c_cell_type_tx;
            case 7 % cluster id
                ci = sce.c_cluster_id;
            case 8 % batch id
                ci = sce.c_batch_id;
            otherwise   % other properties
                ttxt = sce.list_cell_attributes{2 * (indx - 8) - 1};
                ci = sce.list_cell_attributes{2 * (indx - 8)};
        end
        
% ------ move to i_showstate        
        if isempty(ci)
            errordlg("Undefined classification");
            return
        end
        sces = sce.s;
        if isempty(h.ZData)
            sces = sce.s(:, 1:2);
        end        
        [ax, bx] = view();
        h = gui.i_gscatter3(sces, ci, 1);
        view(ax, bx);
        title(sce.title);
% -------- move to i_showstate

        if indx == 5
            hc = colorbar;
            hc.Label.String = ttxt;
        else
            colorbar off;
        end
        [c, cL] = grp2idx(ci);
        sce.c = ci;
        guidata(FigureHandle, sce);
        
        if indx == 5
            answer = questdlg('Compare G1/S/G2M ratios between classes?');
            if ~isequal(answer, 'Yes')
                return
            end
            
            listitems = {'Cluster ID', 'Batch ID', ...
                'Cell Type'};
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select statistics', '', ''}, ...
                'SelectionMode', 'single', ...
                'ListString', listitems);
            if tf2 == 1
                switch indx2
                    case 1 % cluster id
                        thisc = sce.c_cluster_id;
                    case 2 % batch id
                        thisc = sce.c_batch_id;
                    case 3 % cell type
                        thisc = sce.c_cell_type_tx;
                end
            else
                return
            end
            if isempty(thisc)
                errordlg("Undefined classification");
                return
            end
            try
                [A, ~, ~, l] = crosstab(sce.c_cell_cycle_tx, thisc);
                B = A ./ sum(A);
                figure;
                bar(B', 'stacked');
                ylabel(sprintf('%s|', string(l(1:3, 1))));
            catch ME
                rethrow(ME);
            end
        end
    end

    function i_showstate(ci)
        if isempty(ci)
            errordlg("Undefined classification");
            return
        end
        sces = sce.s;
        if isempty(h.ZData)
            sces = sce.s(:, 1:2);
        end        
        [ax, bx] = view();
        h = gui.i_gscatter3(sces, ci, 1);
        view(ax, bx);
        title(sce.title);        
    end

    function DeleteSelectedCells(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return
        end
        answer = questdlg('Delete cells?', '', ...
            'Selected', 'Unselected', 'Cancel', 'Selected');
        if strcmp(answer, 'Cancel')
            return
        end
        if strcmp(answer, 'Unselected')
            ptsSelected = ~ptsSelected;
        end
        answer2 = questdlg(sprintf('Delete %s cells?', ...
            lower(answer)));
        if ~strcmp(answer2, 'Yes')
            return
        end
        i_deletecells(ptsSelected);
        guidata(FigureHandle, sce);
    end

    function i_deletecells(ptsSelected)
        sce = sce.removecells(ptsSelected);
        [c, cL] = grp2idx(sce.c);
        [ax, bx] = view();
        h = gui.i_gscatter3(sce.s, c);
        title(sce.title);
        view(ax, bx);
    end

    function DrawTrajectory(~, ~)
        answer = questdlg('Which method?', 'Select Algorithm', ...
            'splinefit (fast)', 'princurve (slow)', ...
            'splinefit (fast)');
        if strcmp(answer, 'splinefit (fast)')
            dim = 1;
            [t, xyz1] = i_pseudotime_by_splinefit(sce.s, dim, false);
        else
            [t, xyz1] = i_pseudotime_by_princurve(sce.s, false);
        end
        hold on;
        if size(xyz1, 2) >= 3
            plot3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-r', 'linewidth', 2);
            text(xyz1(1, 1), xyz1(1, 2), xyz1(1, 3), 'Start', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            text(xyz1(end, 1), xyz1(end, 2), xyz1(end, 3), 'End', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        elseif size(xyz1, 2) == 2
            plot(xyz1(:, 1), xyz1(:, 2), '-r', 'linewidth', 2);
            text(xyz1(1, 1), xyz1(1, 2), 'Start', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            text(xyz1(end, 1), xyz1(end, 2), 'End', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        end
        hold off;
        
        answerx = questdlg('Save/Update pseudotime T in SCE', ...
            'Save Pseudotime', ...
            'Yes', 'No', 'Yes');
        switch answerx
            case 'Yes'
                tag = sprintf('%s pseudotime', answer);
                % iscellstr(sce.list_cell_attributes(1:2:end))
                i = find(contains(sce.list_cell_attributes(1:2:end), tag));
                if ~isempty(i)
                    sce.list_cell_attributes{i + 1} = t;
                    fprintf('%s is updated.\n', tag);
                else
                    sce.list_cell_attributes{end + 1} = tag;
                    sce.list_cell_attributes{end + 1} = t;
                    fprintf('%s is saved.\n', tag);
                end
                guidata(FigureHandle, sce);
        end
        answer = questdlg('View expression of selected genes', ...
            'Pseudotime Function', ...
            'Yes', 'No', 'Yes');
        switch answer
            case 'Yes'
                r = corr(t, sce.X.', 'type', 'spearman'); % Calculate linear correlation between gene expression profile and T
                [~, idxp] = maxk(r, 4);  % Select top 4 positively correlated genes
                [~, idxn] = mink(r, 3);  % Select top 3 negatively correlated genes
                selectedg = sce.g([idxp idxn]);
                figure;
                i_plot_pseudotimeseries(log2(sce.X + 1), ...
                    sce.g, t, selectedg);
            case 'No'
                return
        end
        
    end

    function RunTrajectoryAnalysis(~, ~)
        answer = questdlg('Run pseudotime analysis (Monocle)?');
        if ~strcmp(answer, 'Yes')
            return
        end
        
        fw = gui.gui_waitbar;
        [t_mono, s_mono] = run.monocle(sce.X);
        gui.gui_waitbar(fw);
        
        answer = questdlg('View Monocle DDRTree?', ...
            'Pseudotime View', ...
            'Yes', 'No', 'Yes');
        switch answer
            case 'Yes'
                [ax, bx] = view();
                cla(hAx);
                sce.s = s_mono;
                sce.c = t_mono;
                [c, cL] = grp2idx(sce.c);
                h = gui.i_gscatter3(sce.s, c);
                title(sce.title);
                view(ax, bx);
                hc = colorbar;
                hc.Label.String = 'Pseudotime';
        end
        
        labels = {'Save pseudotime T to variable named:', ...
            'Save S to variable named:'};
        vars = {'t_mono', 's_mono'};
        values = {t_mono, s_mono};
        export2wsdlg(labels, vars, values);
    end

    function ClusterCellsS(src, ~)
        answer = questdlg('Cluster cells?');
        if ~strcmp(answer, 'Yes')
            return
        end
        
        answer = questdlg('Which method?', 'Select Algorithm', ...
            'kmeans', 'snndpc', 'kmeans');
        if strcmpi(answer, 'kmeans')
            methodtag = "kmeans";
        elseif strcmpi(answer, 'snndpc')
            methodtag = "snndpc";
        else
            return
        end
        i_reclustercells(src, methodtag);
        guidata(FigureHandle, sce);
    end

    function k = i_inputk
        prompt = {'Enter number of clusters K=(2..50):'};
        dlgtitle = 'Input K';
        dims = [1 45];
        definput = {'10'};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer)
            k=[];
            return
        end
        k = round(str2double(cell2mat(answer)));
    end

    function ClusterCellsX(src, ~)
        answer = questdlg('Cluster cells using X?');
        if ~strcmp(answer, 'Yes')
            return
        end
        methodtagv = {'sc3', 'simlr', 'soptsc', 'sinnlrr', 'specter'};
        [indx, tf] = listdlg('PromptString', {'Select clustering program', ...
            '', ''}, 'SelectionMode', 'single', ...
            'ListString', methodtagv);
        if tf == 1
            methodtag = methodtagv{indx};
        else
            return
        end
        i_reclustercells(src, methodtag);
    end

    function i_reclustercells(src, methodtag)
        methodtag = lower(methodtag);
        usingold = false;
        if ~isempty(sce.struct_cell_clusterings.(methodtag))
            answer1 = questdlg(sprintf('Using existing %s clustering?', upper(methodtag)), ...
                '', ...
                'Yes, use existing', 'No, re-compute', 'Cancel', 'Yes, use existing');
            switch answer1
                case 'Yes, use existing'
                    sce.c_cluster_id = sce.struct_cell_clusterings.(methodtag);
                    usingold = true;
                case 'No, re-compute'
                    usingold = false;
                case 'Cancel'
                    return
            end
        end
        if ~usingold
            k = i_inputk;
            if isempty(k), return; end
            if isnan(k) || k < 2 || k > 50
                uiwait(errordlg('Invalid K'));
                return
            end
            fw = gui.gui_waitbar;
            try
                % [sce.c_cluster_id]=sc_cluster_x(sce.X,k,'type',methodtag);
                sce = sce.clustercells(k, methodtag, true);
            catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
                return
            end
            gui.gui_waitbar(fw);
        end
        [c, cL] = grp2idx(sce.c_cluster_id);
        sce.c = c;
        RefreshAll(src, [], true, false);
        guidata(FigureHandle, sce);
    end

    function LabelClusters(src, ~)
        state = src.State;
        if strcmp(state, 'off')
            dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);
        else
            listitems = {'Current Class Type'};
            if ~isempty(sce.c_cluster_id)
                listitems = [listitems, 'Cluster ID'];
            end
            if ~isempty(sce.c_cell_type_tx)
                listitems = [listitems, 'Cell Type'];
            end
            if ~isempty(sce.c_cell_cycle_tx)
                listitems = [listitems, 'Cell Cycle Phase'];
            end
            if ~isempty(sce.c_batch_id)
                listitems = [listitems, 'Batch ID'];
            end
            [indx, tf] = listdlg('PromptString', {'Select statistics', ...
                '', ''}, 'SelectionMode', 'single', 'ListString', listitems);
            if tf ~= 1
                set(src, 'State', 'off');
                return
            end
            switch listitems{indx}
                case 'Cluster ID'
                    cc = sce.c_cluster_id;
                case 'Cell Type'
                    cc = sce.c_cell_type_tx;
                case 'Cell Cycle Phase'
                    cc = sce.c_cell_cycle_tx;
                case 'Batch ID'
                    cc = sce.c_batch_id;
                otherwise
                    cc = [];
            end
            if ~isempty(cc)
                [c, cL] = grp2idx(cc);
                sce.c = c;
            end
            RefreshAll(src, 1, true, false);
            if i_labelclusters
                set(src, 'State', 'on');
            else
                set(src, 'State', 'off');
            end
            guidata(FigureHandle, sce);
            % colormap(lines(min([256 numel(unique(sce.c))])));
        end
    end

    function ShowClustersPop(src, ~)
        answer = questdlg('Show clusters in new figures?');
        if ~strcmp(answer, 'Yes')
            return
        end
        
        cmv = 1:max(c);
        idxx = cmv;
        [cmx] = countmember(cmv, c);
        answer = questdlg('Sort by size of cell groups?');
        if strcmpi(answer, 'Yes')
            [~, idxx] = sort(cmx, 'descend');
        end
        sces = sce.s;
        if isempty(h.ZData)
            sces = sce.s(:, 1:2);
        end
        
        [para] = i_getoldsettings(src);
        figure;
        for k = 1:9
            if k <= max(c)
                subplot(3, 3, k);
                gui.i_gscatter3(sces, c, 3, cmv(idxx(k)));
                title(sprintf('%s\n%d cells (%.2f%%)', ...
                    cL{idxx(k)}, cmx(idxx(k)), ...
                    100 * cmx(idxx(k)) / length(c)));
            end
            colormap(para.oldColorMap);
        end
        
        if ceil(max(c) / 9) == 2
            figure;
            for k = 1:9
                kk = k + 9;
                if kk <= max(c)
                    subplot(3, 3, k);
                    gui.i_gscatter3(sces, c, 3, cmv(idxx(kk)));
                    title(sprintf('%s\n%d cells (%.2f%%)', ...
                        cL{idxx(kk)}, cmx(idxx(kk)), ...
                        100 * cmx(idxx(kk)) / length(c)));
                end
            end
            colormap(para.oldColorMap);
        end
        if ceil(max(c) / 9) > 2
            warndlg('Group(s) #18 and above are not displayed');
        end
    end

    function [txt] = i_myupdatefcnx(~, event_obj)
        % pos = event_obj.Position;
        idx = event_obj.DataIndex;
        txt = cL(c(idx));
    end

    function [isdone] = i_labelclusters(notasking)
        if nargin < 1
            notasking = false;
        end
        isdone = false;
        if ~isempty(cL)
            if notasking
                stxtyes = c;
            else
                answer = questdlg(sprintf('Label %d groups with index or text?', numel(cL)), ...
                    'Select Format', 'Index', 'Text', 'Cancel', 'Text');
                switch answer
                    case 'Text'
                        stxtyes = cL(c);
                    case 'Index'
                        stxtyes = c;
                    otherwise
                        return
                end
            end
            dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);
            
            row = dataTipTextRow('', stxtyes);
            h.DataTipTemplate.DataTipRows = row;
            % h.DataTipTemplate.FontSize = 5;
            for i = 1:max(c)
                idx = find(c == i);
                siv = sce.s(idx, :);
                si = mean(siv, 1);
                [k] = dsearchn(siv, si);
                datatip(h, 'DataIndex', idx(k));
            end
            isdone = true;
        end
    end

    function [para] = i_getoldsettings(src)
        ah = findobj(src.Parent.Parent, 'type', 'Axes');
        ha = findobj(ah.Children, 'type', 'Scatter');
        ha1 = ha(1);
        oldMarker = ha1.Marker;
        oldSizeData = ha1.SizeData;
        oldColorMap = colormap;
        para.oldMarker = oldMarker;
        para.oldSizeData = oldSizeData;
        para.oldColorMap = oldColorMap;
    end

end
