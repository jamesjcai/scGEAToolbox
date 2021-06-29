%function add_toolbar(FigureHandle)
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
    '..','resources', 'list.gif'));
ptImage = ind2rgb(img, map);
pt3.CData = ptImage;
pt3.Tooltip = 'Select a gene to show expression';
pt3.ClickedCallback = @callback_ShowGeneExpr;

pt3a = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'list2.gif'));
ptImage = ind2rgb(img, map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Show cell states';
pt3a.ClickedCallback = @ShowCellStats;

pt3a = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-pointfig.gif'));
ptImage = ind2rgb(img, map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Select cells by class';
pt3a.ClickedCallback = @callback_SelectCellsByClass;

pt3a = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-effects.gif'));
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
    '..','resources', 'plotpicker-glyplot-face.gif'));
ptImage = ind2rgb(img, map);
ptaddcluster.CData = ptImage;
ptaddcluster.Tooltip = 'Add brushed cells to a new cluster';
ptaddcluster.ClickedCallback = @Brushed2NewCluster;

ptmergecluster = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-pzmap.gif'));
ptImage = ind2rgb(img, map);
ptmergecluster.CData = ptImage;
ptmergecluster.Tooltip = 'Merge brushed cells to same cluster';
ptmergecluster.ClickedCallback = @Brushed2MergeClusters;

ptShowClu = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-geoscatter.gif'));
ptImage = ind2rgb(img, map);
ptShowClu.CData = ptImage;
ptShowClu.Tooltip = 'Show clusters individually';
ptShowClu.ClickedCallback = @gui.callback_ShowClustersPop;

ptcluster = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-dendrogram.gif'));
ptImage = ind2rgb(img, map);
ptcluster.CData = ptImage;
ptcluster.Tooltip = 'Clustering using embedding S';
ptcluster.ClickedCallback = @ClusterCellsS;

ptcluster = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-gscatter.gif'));
ptImage = ind2rgb(img, map);
ptcluster.CData = ptImage;
ptcluster.Tooltip = 'Clustering using expression matrix X';
ptcluster.ClickedCallback = @ClusterCellsX;

% -------------

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'brush.gif'));
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
    '..','resources', 'plotpicker-scatterhist.gif'));
ptImage = ind2rgb(img, map);
pt4.CData = ptImage;
pt4.Tooltip = 'Rename cell type';
pt4.ClickedCallback = @RenameCellType;

pt4 = uipushtool(UitoolbarHandle, 'Separator', 'off');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-kagi.gif'));
ptImage = ind2rgb(img, map);
pt4.CData = ptImage;
pt4.Tooltip = 'Marker genes of brushed cells';
pt4.ClickedCallback = @callback_Brush4Markers;

pt4mrkheat = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-plotmatrix.gif'));
ptImage = ind2rgb(img, map);
pt4mrkheat.CData = ptImage;
pt4mrkheat.Tooltip = 'Marker gene heatmap';
pt4mrkheat.ClickedCallback = @callback_MarkerGeneHeatmap;

% --------------------------



ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'IMG00107.GIF'));    % white space
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;



ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'IMG00074.GIF'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Setup R environment';
ptpseudotime.ClickedCallback = @gui.i_setrenv;

ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'IMG00067.GIF'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Run Seurat/R Workflow (R required)';
ptpseudotime.ClickedCallback = @RunSeuratWorkflow;



ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-candle.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Compare Differentiation Potency';
ptpseudotime.ClickedCallback = @callback_ComparePotency;

ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-arxtimeseries.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Run pseudotime analysis (Monocle)';
ptpseudotime.ClickedCallback = @RunTrajectoryAnalysis;

ptpseudotime = uipushtool(defaultToolbar, ...
    'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-comet.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Plot pseudotime trajectory';
ptpseudotime.ClickedCallback = @DrawTrajectory;

ptpseudotime = uipushtool(defaultToolbar, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-priceandvol.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Compare Gene Expression between Classes';
ptpseudotime.ClickedCallback = @callback_CompareGeneBtwCls;

pt4 = uipushtool(defaultToolbar, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-boxplot.gif'));
ptImage = ind2rgb(img, map);
pt4.CData = ptImage;
pt4.Tooltip = 'Compare 2 groups (DE analysis)';
pt4.ClickedCallback = @callback_DEGene2Groups;

ptpseudotime = uipushtool(defaultToolbar, ...
    'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-andrewsplot.gif'));
ptImage = ind2rgb(img, map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Function enrichment of HVG genes';
ptpseudotime.ClickedCallback = @callback_GSEA_HVGs;

ptnetwork = uipushtool(defaultToolbar, ...
    'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'noun_Network_691907.gif'));
ptImage = ind2rgb(img, map);
ptnetwork.CData = ptImage;
ptnetwork.Tooltip = 'Build gene regulatory network';
ptnetwork.ClickedCallback = @callback_BuildGeneNetwork;

ptnetwork = uipushtool(defaultToolbar, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'noun_Deep_Learning_2424485.gif'));
ptImage = ind2rgb(img, map);
ptnetwork.CData = ptImage;
ptnetwork.Tooltip = 'Compare two scGRNs';
ptnetwork.ClickedCallback = @callback_CompareGeneNetwork;

pt2 = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-qqplot.gif'));
ptImage = ind2rgb(img, map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delete selected cells';
pt2.ClickedCallback = @DeleteSelectedCells;

pt = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, '..','resources', 'export.gif'));
ptImage = ind2rgb(img, map);
pt.CData = ptImage;
pt.Tooltip = 'Export & save data';
pt.ClickedCallback = @callback_SaveX;

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-geobubble.gif'));
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Embedding';
pt5.ClickedCallback = @EmbeddingAgain;

% run(fullfile(mfolder,'+gui','add_toolbar.m'))
% pt5 = uipushtool(UitoolbarHandle, 'Separator', 'off');
% [img, map] = imread(fullfile(mfolder, ...
%     '..','resources', 'multiscale.gif'));
% ptImage = ind2rgb(img, map);
% pt5.CData = ptImage;
% pt5.Tooltip = 'Run Seurat/R Workflow (R required)';
% pt5.ClickedCallback = @RunSeuratWorkflow;

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-image.gif'));      % plotpicker-pie
% map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;     % Convert white pixels => transparent background
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Switch 2D/3D';
pt5.ClickedCallback = @Switch2D3D;

pt5pickmk = uipushtool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-rose.gif'));  % plotpicker-pie
ptImage = ind2rgb(img, map);
pt5pickmk.CData = ptImage;
pt5pickmk.Tooltip = 'Switch scatter plot marker type';
pt5pickmk.ClickedCallback = @callback_PickPlotMarker;

pt5pickcl = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-compass.gif'));  % plotpicker-pie
ptImage = ind2rgb(img, map);
pt5pickcl.CData = ptImage;
pt5pickcl.Tooltip = 'Switch color maps';
pt5pickcl.ClickedCallback = {@gui.callback_PickColorMap, ...
    10};
%    numel(unique(c))};

pt5 = uipushtool(UitoolbarHandle, 'Separator', 'off');
[img, map] = imread(fullfile(mfolder, ...
    '..','resources', 'plotpicker-geobubble2.gif'));
ptImage = ind2rgb(img, map);
pt5.CData = ptImage;
pt5.Tooltip = 'Refresh';
pt5.ClickedCallback = @RefreshAll;

gui.add_3dcamera(defaultToolbar, 'AllCells');