function varargout = scgeatoolx(sce, varargin)

if usejava('jvm') && ~feature('ShowFigureWindows')
    error('MATLAB is in a text mode. This function requires a GUI-mode.');
end

if nargin < 1
   sce = SingleCellExperiment; 
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
addOptional(p, 'callinghandle', []);
parse(p, sce, varargin{:});
callinghandle = p.Results.callinghandle;


c_in = p.Results.c;
s_in = p.Results.s;
methodid = p.Results.methodid;
f_traj = [];   % trajectory curve

ax = [];
bx = [];
pushbuttonV = [];

tmpcelltypev = cell(sce.NumCells, 1);

if ~isempty(c_in), sce.c = c_in; end
if ~isempty(s_in), sce.s = s_in; end

[c, cL] = grp2idx(sce.c);

if ~(ismcc || isdeployed)
    tagx = 'on';
else
    tagx = 'off';
end


FigureHandle = figure('Name', 'SCGEATOOL :: Single-Cell Gene Expression Analysis Tool', ...
    'position', round(1.25*[0, 0, 560, 420]), ...
    'visible', 'off', 'NumberTitle', tagx);
movegui(FigureHandle, 'center');

% b=uipanel(FigureHandle,'Title','B','BackgroundColor','cyan');
% b.Position = [0.18 0.40 0.30 0.35];

% set(findall(FigureHandle, 'ToolTipString', 'Insert Colorbar'), 'Visible', 'Off')
% set(findall(FigureHandle, 'ToolTipString', 'Insert Legend'), 'Visible', 'Off')
% set(findall(FigureHandle, 'ToolTipString', 'Print Figure'), 'Visible', 'Off')

pushbuttonV = findall(FigureHandle, 'ToolTipString', 'Insert Colorbar');
set(pushbuttonV,'Separator','off');
%pushbuttonV = [pushbuttonV; findall(FigureHandle, 'ToolTipString', 'Insert Legend')];
%pushbuttonV = [pushbuttonV; findall(FigureHandle, 'ToolTipString', 'Print Figure')];
delete(findall(FigureHandle, 'ToolTipString', 'Print Figure'));
delete(findall(FigureHandle, 'ToolTipString', 'Insert Legend'));
delete(findall(FigureHandle, 'ToolTipString', 'Link/Unlink Plot'));
delete(findall(FigureHandle, 'ToolTipString', 'Edit Plot'));
delete(findall(FigureHandle, 'ToolTipString', 'Open Property Inspector'));
delete(findall(FigureHandle, 'ToolTipString', 'Save Figure'));
delete(findall(FigureHandle, 'ToolTipString', 'New Figure'));
delete(findall(FigureHandle, 'ToolTipString', 'Open File'));



%a=findall(FigureHandle,'ToolTipString','New Figure');
%a.ClickedCallback = @__;
%a=findall(f,'tag','figMenuFile');

% https://undocumentedmatlab.com/articles/customizing-standard-figure-toolbar-menubar

% delete(findall(FigureHandle,'tag','figMenuEditClearWorkspace'));
% delete(findall(FigureHandle,'tag','figMenuEditClearCmdHistory'));
% delete(findall(FigureHandle,'tag','figMenuEditClearCmdWindow'));
% delete(findall(FigureHandle,'tag','figMenuEditClearFigure'));
% delete(findall(FigureHandle,'tag','figMenuEditFindFiles'));
% delete(findall(FigureHandle,'tag','figMenuEditSelectAll'));
% delete(findall(FigureHandle,'tag','figMenuEditGCO'));
% delete(findall(FigureHandle,'tag','figMenuEditGCA'));
% delete(findall(FigureHandle,'tag','figMenuEditGCF'));
delete(allchild(findall(FigureHandle, 'tag', 'figMenuView')))

delete(findall(FigureHandle,'tag','figMenuInsertAxes'));
delete(findall(FigureHandle,'tag','figMenuInsertLight'));
delete(findall(FigureHandle, 'tag', 'figMenuWindow'));
delete(findall(FigureHandle, 'tag', 'figMenuDesktop'));
delete(findall(FigureHandle, 'tag', 'figMenuUpdateFileNew'));
set(findall(FigureHandle, 'tag', 'figMenuOpen'),'MenuSelectedFcn', @in_sc_openscedlg);
delete(findall(FigureHandle, 'tag', 'figMenuFileSaveAs'));
delete(findall(FigureHandle, 'tag', 'figMenuFileSaveWorkspaceAs'));
delete(findall(FigureHandle, 'tag', 'figMenuFilePreferences'));
delete(findall(FigureHandle, 'tag', 'figMenuFileExportSetup'));
set(findall(FigureHandle, 'tag', 'figMenuFileSave'),'Text','&Save SCE...','MenuSelectedFcn', @callback_SaveX);
delete(findall(FigureHandle, 'tag', 'figMenuGenerateCode'));
delete(findall(FigureHandle, 'tag', 'figMenuFileImportData'));

% set(findall(FigureHandle, 'tag', 'figMenuFileImportData'),'MenuSelectedFcn', @in_GEOAccessionToSCE,...
%     'Text','Import Data Using GEO Accession...','Separator','on');
set(findall(FigureHandle,'tag','figMenuFilePrintPreview'),'Separator','on');

m_file=findall(FigureHandle,'tag','figMenuFile');
in_addmenu(m_file, 1, @in_SimulateSCE, 'Simulate Data [PMID:27122128]...');
in_addmenu(m_file, 1, {@gui.i_savemainfig, 3}, 'Save Figure to PowerPoint File...');
in_addmenu(m_file, 0, {@gui.i_savemainfig, 2}, 'Save Figure as Graphic File...');
in_addmenu(m_file, 0, {@gui.i_savemainfig, 1}, 'Save Figure as SVG File...');


m_edit=findall(FigureHandle,'tag','figMenuEdit');
delete(allchild(m_edit));
in_addmenu(m_edit, 0, @in_SelectCellsByQC, 'Filter genes and cells...');
in_addmenu(m_edit, 0, @in_Brushed2NewCluster, 'Add brushed cells to a new group');
in_addmenu(m_edit, 0, @in_Brushed2MergeClusters, 'Merge brushed cells to same group');
in_addmenu(m_edit, 0, @in_RenameCellTypeBatchID, 'Rename Cell Type or Batch ID...');
in_addmenu(m_edit, 1, @gui.callback_SplitAtacGex, 'Split Multiome ATAC+GEX Matrix...');
in_addmenu(m_edit, 1, {@in_MergeSCEs, 1}, 'Merge SCE Variables in Workspace...');
in_addmenu(m_edit, 0, {@in_MergeSCEs, 2}, 'Merge SCE Data Files...');
in_addmenu(m_edit, 1, @in_AddEditCellAttribs, 'Add/Edit Cell Attributes...');
in_addmenu(m_edit, 0, @in_ExportCellAttribTable, 'Export Cell Attribute Table...');
in_addmenu(m_edit, 1, {@in_MergeCellSubtypes, 1}, 'Import Cell Annotation from SCE in Workspace...');
in_addmenu(m_edit, 0, {@in_MergeCellSubtypes, 2}, 'Import Cell Annotation from SCE Data File...');
in_addmenu(m_edit, 1, @gui.callback_SelectCellsByMarker, 'Extract Cells by Marker (+/-) Expression...');
in_addmenu(m_edit, 0, @in_MergeSubCellTypes, 'Merge Subclusters of the Same Cell Type');
%i_addmenu(m_edit,0,@callback_TwoGeneCooccurrenceTest,'Two-Gene Cooccurrence Test...');
%i_addmenu(m_edit,0,@AnnotateSubGroup,'Annotate Cell Subgroups...');
in_addmenu(m_edit, 1, @in_WorkonSelectedGenes, 'Select Top n  Highly Variable Genes (HVGs) to Work on...');
in_addmenu(m_edit, 0, @in_SubsampleCells, 'Subsample 50% Cells to Work on...');
in_addmenu(m_edit, 1, @gui.callback_SelectCellsByClass, 'Select Cells...');

%i_addmenu(m_exp,0,{@MergeCellSubtypes,1,true},'Import All Cell Annotation from SCE in Workspace...');
%i_addmenu(m_exp,0,{@MergeCellSubtypes,2,true},'Import All Cell Annotation from SCE Data File...');

m_view=findall(FigureHandle,'tag','figMenuView');
in_addmenu(m_view, 0, @gui.callback_ShowGeneExpr, 'Gene Expression...');
in_addmenu(m_view, 0, @in_ShowCellStates, 'Cell States...');
in_addmenu(m_view, 0, @in_labelcellgroups, 'Cell Groups...');
in_addmenu(m_view, 0, @gui.callback_MultiGroupingViewer, 'Multi-Grouping View...');
in_addmenu(m_view, 0, @gui.callback_CrossTabulation, 'Cross Tabulation');
in_addmenu(m_view, 1, @gui.callback_ViewMetaData, 'View Metadata...');
in_addmenu(m_view, 1, @gui.callback_ShowHgBGeneExpression, 'Hemoglobin (Hgb) Genes Expression...');
in_addmenu(m_view, 0, @gui.callback_ShowMtGeneExpression, 'Mitochondrial (Mt-) Genes Expression...');
in_addmenu(m_view, 1, @in_DrawKNNNetwork, 'Cell kNN Network...');
in_addmenu(m_view, 0, @in_DrawTrajectory, 'Cell Trajectory...');
in_addmenu(m_view, 1, @gui.callback_CloseAllOthers, 'Close All Other Figures');
in_addmenu(m_view, 0, @in_RefreshAll, 'Refresh Current View');

m_plot = findall(FigureHandle, 'tag', 'figMenuTools');
set(m_plot,'Text','&Plot')
delete(allchild(m_plot));
set(findall(FigureHandle,'tag','figMenuInsert'),'Parent', m_plot,'Separator','off');
set(findall(FigureHandle,'tag','figMenuEditCopyFigure'),'Parent', m_plot);
set(findall(FigureHandle,'tag','figMenuEditCopyOptions'),'Parent', m_plot);
in_addmenu(m_plot,1,@gui.callback_PickColorMap,'Next Colormap');
set(findall(FigureHandle,'tag','figMenuEditColormap'),'Parent', m_plot,'Text','Colormap Editor...');
in_addmenu(m_plot,0,@gui.callback_PickPlotMarker,'Next Marker Type');
in_addmenu(m_plot,1,@gui.callback_Violinplot,'Gene Violin Plot...');
in_addmenu(m_plot,0,@gui.callback_DrawDotplot,'Gene Dot Plot...');
in_addmenu(m_plot,0,@gui.callback_GeneHeatMap,'Gene Heatmap...');

m_exp = uimenu(FigureHandle, 'Text', '&Tools', 'Accelerator', 'T');
% m_exp2 = uimenu(m_exp,'Text','sc&Tenifold Suite','Accelerator','T');
% i_addmenu(m_exp2,1,@callback_scTenifoldNet1,'scTenifoldNet - GRN Construction 🐢🐢 ...');
% i_addmenu(m_exp2,0,@callback_scTenifoldNet2,'scTenifoldNet - GRN Comparison 🐢🐢🐢 ...');
% i_addmenu(m_exp2,0,@callback_scTenifoldKnk1,'scTenifoldKnk - Virtual KO of a Gene 🐢 ...');
% i_addmenu(m_exp2,0,@callback_scTenifoldXct,'scTenifoldXct - Cell-Cell Interactions 🐢 ...');
%i_addmenu(m_exp,0,@callback_ShowPseudoTimeGenes,'Show Genes with Expression Varies with Pseudotime...');
%i_addmenu(m_exp,0,@callback_DetectCellularCrosstalk,'Ligand-Receptor Mediated Intercellular Crosstalk...');
%i_addmenu(m_exp,0,@AnnotateSubTypes,'Assign Subtypes of Cells (Neurons or T Cells)...');
in_addmenu(m_exp, 0, @in_CompareGeneBtwCls, 'Cell Score Analysis...');
in_addmenu(m_exp, 0, @gui.callback_GetCellSignatureMatrix, 'Cell State Analysis...');

in_addmenu(m_exp, 1, @in_EnrichrHVGs, 'HVG Functional Enrichment Analysis...');
in_addmenu(m_exp, 1, @gui.callback_DEGene2Groups, 'Differential Expression (DE) Analysis...');
in_addmenu(m_exp, 0, @gui.callback_DPGene2Groups, 'Differential Program (DP) Analysis...');
in_addmenu(m_exp, 0, @gui.callback_DEGene2GroupsBatch, 'Differential Expression (DE) Analysis in Cell Type Batch Mode...');
in_addmenu(m_exp, 0, @gui.callback_DPGene2GroupsBatch, 'Differential Program (DP) Analysis in Cell Type Batch Mode...');

in_addmenu(m_exp, 1, @gui.callback_CalculateGeneStats, 'Calculate Gene Expression Statistics...');
in_addmenu(m_exp, 0, @gui.callback_CellCycleLibrarySize, 'Library Size of Cell Cycle Phases...');
in_addmenu(m_exp, 0, @gui.callback_CellCycleAssignment, 'Assign Cell Cycle Phase...');
%i_addmenu(m_exp,0,@gui.callback_TCellExhaustionScores,'T Cell Exhaustion Score...');
in_addmenu(m_exp, 1, {@in_DetermineCellTypeClustersGeneral, false}, 'Annotate Cell Type Using Customized Markers...');
in_addmenu(m_exp, 0, @in_SubtypeAnnotation, 'Annotate Cell Subtype...');

in_addmenu(m_exp, 1, @in_SingleClickSolution, 'Single Click Solution (from Raw Data to Annotation)...');


m_net = uimenu(FigureHandle, 'Text', '&Network', 'Accelerator', 'N');
in_addmenu(m_net, 0, @in_Select5000Genes, 'Remove Less Informative Genes to Reduce Gene Space...');
in_addmenu(m_net, 0, @gui.i_setnetwd, 'Set Network Analysis Working Root Directory...');
in_addmenu(m_net, 1, @gui.callback_BuildGeneNetwork, 'Build GRN with Selected Genes...');
in_addmenu(m_net, 0, @gui.callback_CompareGeneNetwork, 'Build & Compare Two GRNs...');
in_addmenu(m_net, 1, {@in_scTenifoldNet,1}, 'Construct GRN using PC Regression [PMID:33336197] 🐢...');
%in_addmenu(m_net, 1, @callback_scPCNet1, 'GRN Construction - PC Regression (w/o subsampling) [PMID:33336197] 🐢...');
%in_addmenu(m_net, 0, @callback_scTenifoldNet1, 'GRN Construction - PC Regression (w/ subsampling) [PMID:33336197] 🐢🐢 ...');
in_addmenu(m_net, 0, {@in_scTenifoldNet,2}, 'Construct & Compare GRNs (scTenifoldNet Analysis) [PMID:33336197] 🐢...');

%in_addmenu(m_net, 1, @callback_scTenifoldNet2lite, 'GRN Comparison - scTenifoldNet (w/o subsampling) [PMID:33336197] 🐢🐢 ...');
%in_addmenu(m_net, 0, @callback_scTenifoldNet2, 'GRN Comparison - scTenifoldNet (w/ subsampling) [PMID:33336197] 🐢🐢🐢 ...');

in_addmenu(m_net, 1, @gui.callback_scTenifoldKnk1, 'Virtual Gene KO - scTenifoldKnk [PMID:35510185] 🐢🐢 ...');
in_addmenu(m_net, 0, @gui.callback_VirtualKOGenKI, 'Virtual Gene KO - GenKI [PMID:37246643] (Python Required) 🐢🐢 ...');
in_addmenu(m_net, 1, @gui.callback_scTenifoldXct, 'Cell-Cell Interactions (CCIs) - scTenifoldXct [PMID:36787742] 🐢🐢 ...');
in_addmenu(m_net, 0, @gui.callback_scTenifoldXct2, 'Differential CCIs - scTenifoldXct [PMID:36787742] 🐢🐢🐢 ...');

m_ext = uimenu(FigureHandle, 'Text', 'E&xternal', 'Accelerator', 'x');
in_addmenu(m_ext, 0, @gui.i_setrenv, 'Check R Environment');
in_addmenu(m_ext, 0, @gui.i_setpyenv, 'Check Python Environment');
in_addmenu(m_ext, 0, @gui.i_setextwd, 'Set External Program Working Root Directory...');

in_addmenu(m_ext, 1, @in_DecontX, 'Detect Ambient RNA Contamination (DecontX/R) [PMID:32138770]...');
%i_addmenu(m_ext,0,@callback_SingleRCellType,'SingleR Cell Type Annotation (SingleR/R required)...');
%i_addmenu(m_ext,0,@callback_RevelioCellCycle,'Revelio Cell Cycle Analysis (Revelio/R required)...');
% i_addmenu(m_ext,0,@callback_RunSeuratSCTransform,'Run Seurat/R SCTransform (Seurat/R required)...');
in_addmenu(m_ext, 0, @in_RunSeuratWorkflow, 'Run Seurat/R Workflow (Seurat/R) [PMID:25867923]...');
%i_addmenu(m_ext,0,@callback_RunMonocle2,'Pseudotime Analysis (Monocle2/R) [PMID:28825705]...');
in_addmenu(m_ext, 0, @gui.callback_RunMonocle3, 'Pseudotime Analysis (Monocle3/R) [PMID:28825705]...');
in_addmenu(m_ext, 1, @gui.callback_MELDPerturbationScore, 'MELD Perturbation Score (MELD/Py) [PMID:33558698]...');
in_addmenu(m_ext, 0, {@in_SubsampleCells, 2}, 'Geometric Sketching (geosketch/Py) [PMID:31176620]...');
in_addmenu(m_ext, 0, @in_HarmonyPy, 'Batch Integration (Harmony/Py) [PMID:31740819]...');
in_addmenu(m_ext, 0, @in_DoubletDetection, 'Detect Doublets (Scrublet/Py) [PMID:30954476]...');
in_addmenu(m_ext, 0, @in_RunDataMapPlot, 'Run DataMapPlot (datamapplot/Py)...');
in_addmenu(m_ext, 1, @gui.callback_ExploreCellularCrosstalk, 'Talklr Intercellular Crosstalk [DOI:10.1101/2020.02.01.930602]...');

% in_addmenu(m_ext, 0, @gui.callback_CompareGCLBtwCls, 'Differential GCL Analysis [PMID:33139959]🐢🐢 ...');
% in_addmenu(m_ext, 0, @gui.callback_DiffTFActivity, 'Differential TF Activity Analysis...');

delete(findall(FigureHandle, 'tag', 'figMenuHelp'));
m_help = uimenu(FigureHandle, 'Text', '&Help', 'Accelerator', 'H');
in_addmenu(m_help, 0, {@(~, ~) web('https://scgeatoolbox.readthedocs.io/en/latest/')}, 'Online Documentation...');
in_addmenu(m_help, 0, {@(~, ~) web('https://scgeatoolbox.readthedocs.io/en/latest/quick_installation.html')}, 'Quick Installation...');
in_addmenu(m_help, 1, {@(~, ~) web('https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox-single-cell-gene-expression-analysis-toolbox')}, 'View scGEAToolbox on File Exchange...');
in_addmenu(m_help, 0, {@(~, ~) web('https://scgeatool.github.io/')}, 'Visit SCGEATOOL-Standalone Website...');
in_addmenu(m_help, 1, {@(~, ~) web('https://pubmed.ncbi.nlm.nih.gov/31697351/')}, 'Cite scGEAToolbox Paper...');
in_addmenu(m_help, 0, {@(~, ~) web('https://scholar.google.com/scholar?cites=4661048952867744439&as_sdt=5,44&sciodt=0,44&hl=en')}, 'Show Papers Citing scGEAToolbox...');
in_addmenu(m_help, 1, {@(~, ~) web('https://matlab.mathworks.com/open/github/v1?repo=jamesjcai/scGEAToolbox&file=online_landing.m')}, 'Open in MATLAB Online...');
in_addmenu(m_help, 1, @callback_CheckUpdates, 'Check for Updates...');
in_addmenu(m_help, 0, {@(~, ~) web('https://scgeatoolbox.readthedocs.io/en/latest/license.html')}, 'License Agreement');

hAx = axes('Parent', FigureHandle,'Visible','off');
if ~isempty(sce) && sce.NumCells>0
    h = gui.i_gscatter3(sce.s, c, methodid, 1, hAx);
else
    h = [];
end
title(hAx, sce.title);
subtitle(hAx,'[genes x cells]');

dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

mfolder = fileparts(mfilename('fullpath'));

DeftToolbarHandle = findall(FigureHandle, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
MainToolbarHandle = uitoolbar('Parent', FigureHandle);
set(MainToolbarHandle, 'Tag', 'MainToolBar', 'HandleVisibility', 'off', 'Visible', 'off');
UserToolbarHandle = uitoolbar('Parent', FigureHandle);
set(UserToolbarHandle, 'Tag', 'UserToolBar', 'HandleVisibility', 'off', 'Visible', 'off');

% i_addbutton(1,1,@callback_Brush4Markers,"icon-mat-filter-1-20.gif","Marker genes of brushed cells");
% i_addbutton_toggle(1,0,{@togglebtfun,@turnoffuserguiding,"icon-mat-blur-off-10.gif", ...
%    "icon-mat-blur-on-10.gif",false},"Turn on/off user onboarding toolbar");

in_addbuttontoggle(1, 0, {@in_togglebtfun, @in_turnoffuserguiding, ...
    "icon-mat-unfold-more-10.gif", ...
    "icon-mat-unfold-less-10.gif", false, ...
    "Turn on/off user onboarding toolbar"});
in_addbuttonpush(1, 0, @gui.callback_ShowGeneExpr, "list.gif", "Select genes to show expression")
in_addbuttonpush(1, 0, @in_ShowCellStates, "list2.gif", "Show cell state")
in_addbuttonpush(1, 0, @in_SelectCellsByQC, "plotpicker-effects.gif", "Filter genes and cells")

in_addbuttontoggle(1, 1, {@in_togglebtfun, @in_labelcellgroups, ...
    "icon-fa-tag-10b.gif", "icon-fa-tags-10b.gif", ...
    false, "Label cell groups"});
in_addbuttonpush(1, 0, @in_Brushed2NewCluster, "plotpicker-glyplot-face.gif", "Add brushed cells to a new group")
in_addbuttonpush(1, 0, @in_Brushed2MergeClusters, "plotpicker-pzmap.gif", "Merge brushed cells to same group")
in_addbuttonpush(1, 0, @in_RenameCellTypeBatchID, "plotpicker-scatterhist.gif", "Rename cell type or batch ID");
in_addbuttonpush(1, 0, @in_call_scgeatool, "IMG00107.GIF", " ");
in_addbuttonpush(1, 1, @in_ClusterCellsS, "plotpicker-dendrogram.gif", "Clustering using cell embedding (S)")
in_addbuttonpush(1, 0, @in_ClusterCellsX, "icon-mw-cluster-10.gif", "Clustering using expression matrix (X)")
in_addbuttonpush(1, 1, {@in_DetermineCellTypeClustersGeneral, true}, "plotpicker-contour.gif", "Assign cell types to groups")
in_addbuttonpush(1, 0, @in_Brush4Celltypes, "brush.gif", "Assign cell type to selected cells");
% i_addbutton(1,0,@ShowCellStemScatter,"IMG00067.GIF","Stem scatter plot");
in_addbuttonpush(1, 1, @gui.callback_Brush4Markers, "plotpicker-kagi.gif", "Marker genes of brushed cells");
in_addbuttonpush(1, 0, @gui.callback_FindAllMarkers, "plotpicker-plotmatrix.gif", "Marker gene heatmap");
in_addbuttonpush(1, 0, @in_call_scgeatool, "IMG00107.GIF", " ");
in_addbuttonpush(1, 1, @gui.callback_ShowClustersPop, "plotpicker-geoscatter.gif", "Show cell clusters/groups individually");
in_addbuttonpush(1, 0, @gui.callback_SelectCellsByClass, "plotpicker-pointfig.gif", "Select cells by class");
in_addbuttonpush(1, 0, @in_DeleteSelectedCells, "plotpicker-qqplot.gif", "Delete selected cells");
in_addbuttonpush(1, 0, @callback_SaveX, "export.gif", "Export & save data");
in_addbuttonpush(1, 1, @in_EmbeddingAgain, "plotpicker-geobubble.gif", "Embedding (tSNE, UMP, PHATE)");
in_addbuttonpush(1, 0, @in_Switch2D3D, "plotpicker-image.gif", "Switch 2D/3D");
in_addbuttonpush(1, 1, @gui.callback_CloseAllOthers, "icon-fa-cut-10.gif", "Close all other figures");
in_addbuttonpush(1, 0, @gui.callback_PickPlotMarker, "plotpicker-rose.gif", "Switch scatter plot marker type");
in_addbuttonpush(1, 0, @gui.callback_PickColorMap, "plotpicker-compass.gif", "Pick new color map");
in_addbuttonpush(1, 0, @in_RefreshAll, "icon-mat-refresh-20.gif", "Refresh");


%in_addbuttonpush(0, 0, @in_call_scgeatool, "IMG00107.GIF", " ");
%i_addbutton(0,0,@callback_CalculateCellScores,"cellscore2.gif","Calculate cell scores from list of feature genes")
%i_addbutton(0,0,@callback_ComparePotency,"plotpicker-candle.gif","Compare differentiation potency between groups");


in_addbuttonpush(0, 1, @gui.callback_MultiGroupingViewer, "plotpicker-arxtimeseries.gif", "Multi-grouping View...");
in_addbuttonpush(0, 0, @gui.callback_CrossTabulation, "plotpicker-comet.gif", "Cross tabulation");

in_addbuttonpush(0, 1, @gui.callback_Violinplot, "violinplot.gif", "Gene Violin Plot...");
in_addbuttonpush(0, 0, @gui.callback_DrawDotplot, "icon-mat-blur-linear-10.gif", "Gene Dot Plot...");
in_addbuttonpush(0, 0, @gui.callback_GeneHeatMap, "icon-mat-apps-20.gif", "Gene Heatmap...");

in_addbuttonpush(0, 0, @in_call_scgeatool, "IMG00107.GIF", " ");
in_addbuttonpush(0, 1, @in_CompareGeneBtwCls, "cellscore2.gif", "Cell score analysis--obtaining gene signature score for each cell");
in_addbuttonpush(0, 0, @gui.callback_GetCellSignatureMatrix, "icon-fa-connectdevelop-20.gif", "Cell state analysis--obtaining multiple gene signature scores to reveal functional state of cells");
in_addbuttonpush(0, 1, @in_EnrichrHVGs, "plotpicker-andrewsplot.gif", "Functional enrichment analysis with HVGs");
in_addbuttonpush(0, 0, @gui.callback_DEGene2Groups, "plotpicker-boxplot.gif", "Differential expression (DE) analysis)");
in_addbuttonpush(0, 0, @gui.callback_DPGene2Groups, "plotpicker_noisepsd.gif", "Differential program (DP) analysis)");
% fullfile(matlabroot,'toolbox','matlab','icons','HDF_grid.gif')
in_addbuttonpush(0, 1, @gui.callback_BuildGeneNetwork, "noun_Network_691907.gif", "Build gene regulatory network");
in_addbuttonpush(0, 0, @gui.callback_CompareGeneNetwork, "noun_Deep_Learning_2424485.gif", "Compare two scGRNs");
in_addbuttonpush(0, 1, {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');

pushbuttonV=[pushbuttonV; gui.add_3dcamera(DeftToolbarHandle, 'AllCells')];

in_addbuttonpush(2, 0, @in_turnonuserguiding, "icon-fa-thumb-tack-10.gif", ...
    "Turn on user guiding toolbar");

in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_SelectCellsByQC, ...
    "icon-mat-filter-1-10.gif", "plotpicker-effects.gif", ...
    true, "Filter genes and cells"});
in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_EmbeddingAgain, ...
    "icon-mat-filter-2-10.gif", "plotpicker-geobubble.gif", ...
    true, "Embedding (tSNE, UMP, PHATE)"});
in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_ClusterCellsS, ...
    "icon-mat-filter-3-10.gif", "plotpicker-dendrogram.gif", ...
    true, "Clustering using embedding S"});
in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_DetermineCellTypeClustersGeneral, ...
    "icon-mat-filter-4-10.gif", "plotpicker-contour.gif", true, ...
    "Assign cell types to groups"});
in_addbuttontoggle(2, 0, {@in_togglebtfun, @callback_SaveX, ...
    "icon-mat-filter-5-10.gif", "export.gif", ...
    true, "Export & save data"});


% handles = guihandles( FigureHandle );
% guidata( FigureHandle, handles );
if ~isempty(c)
    kc = numel(unique(c));
    colormap(pkg.i_mycolorlines(kc));
end

if ~isempty(sce) && sce.NumCells>0
    in_EnDisableMenu('on');
else
    in_EnDisableMenu('off');
end
drawnow;
if ~isempty(sce) && sce.NumCells>0
    hAx.Visible="on";
end
set(FigureHandle, 'visible', 'on');

in_fixfield('tsne','tsne3d');
in_fixfield('umap','umap3d');
in_fixfield('phate','phate3d');
in_fixfield('metaviz','metaviz3d');

avx = fieldnames(sce.struct_cell_embeddings);
bvx = fieldnames(pkg.e_makeembedstruct);
cvx = setdiff(bvx,avx);
for kx=1:length(cvx)
    sce.struct_cell_embeddings = setfield(sce.struct_cell_embeddings,cvx{kx},[]);
end
sce.struct_cell_embeddings = orderfields(sce.struct_cell_embeddings);

% guidata(FigureHandle, sce);
% setappdata(FigureHandle,'cL',cL);
set(FigureHandle, 'CloseRequestFcn', @in_closeRequest);

if nargout > 0
    varargout{1} = FigureHandle;
end

if ~ispref('scgeatoolbox', 'useronboardingtoolbar')
    gui.gui_userguidingpref(true);
    setpref('scgeatoolbox', 'useronboardingtoolbar', true);
end
showuseronboarding = getpref('scgeatoolbox', 'useronboardingtoolbar');
if ~showuseronboarding
    set(UserToolbarHandle, 'Visible', 'off');
end

if isempty(sce) || sce.NumCells==0
    %[sce]=gui.sc_openscedlg;
    in_sc_openscedlg;
end
%guidata(FigureHandle, sce);
%in_RefreshAll;

% ----------------------------------
% ----------------------------------
    function in_sc_openscedlg(~, ~)
        [sce] = gui.sc_openscedlg;
        if ~isempty(sce) && sce.NumCells > 0
            guidata(FigureHandle, sce);
        end
        in_RefreshAll([], [], false, false);
    end

    function in_fixfield(oldf, newf)
        if ~isfield(sce.struct_cell_embeddings,newf) && isfield(sce.struct_cell_embeddings,oldf)
            if ~isempty(sce.struct_cell_embeddings.(oldf))
                if size(sce.struct_cell_embeddings.(oldf),2) == 3
                    sce.struct_cell_embeddings.(newf) = sce.struct_cell_embeddings.(oldf);
                    sce.struct_cell_embeddings = rmfield(sce.struct_cell_embeddings,oldf);
                end
            end
        end
        if isfield(sce.struct_cell_embeddings, oldf)
            if isempty(sce.struct_cell_embeddings.(oldf))
                sce.struct_cell_embeddings = rmfield(sce.struct_cell_embeddings,oldf);
            end
        end
    end

    function in_CompareGeneBtwCls(src,events)
        gui.callback_CompareGeneBtwCls(src,events);
        sce = guidata(FigureHandle);
    end

    function in_turnonuserguiding(~, ~)
        % setpref('scgeatoolbox','useronboardingtoolbar',true);
        % set(UserToolbarHandle, 'Visible', 'on');
        Button=gui.gui_userguidingpref(false);
        switch Button
            case 'Yes'
            case 'No'
                in_turnoffuserguiding;
            case 'Cancel'
        end
    end

    function in_turnoffuserguiding(~, ~)
        % getpref('scgeatoolbox','useronboardingtoolbar');
        if get(UserToolbarHandle, 'Visible') == "off"
            askpref = true;
        else
            askpref = false;
        end
        if showuseronboarding
            set(UserToolbarHandle, 'Visible', 'off');
        else
            set(UserToolbarHandle, 'Visible', 'on');
        end
        showuseronboarding = ~showuseronboarding;

        if askpref
            %  gui.gui_userguidingpref(false);
            %answer=questdlg('Show User Onboarding Toolbar again next time?','');
            %switch answer
            %    case 'Yes'
            %        setpref('scgeatoolbox','useronboardingtoolbar',true);
            %    case 'No'
            %        setpref('scgeatoolbox','useronboardingtoolbar',false);
            %end
        end
    end

    function in_addmenu(menuHdl, sepTag, callbackFnc, tooltipTxt)
        if ischar(callbackFnc) || isstring(callbackFnc)
            callbackFnc = str2func(callbackFnc);
        end
        if sepTag == 1
            septag = 'on';
        else
            septag = 'off';
        end
        uimenu(menuHdl, 'Text', tooltipTxt, ...
            'Separator', septag, ...
            'Callback', callbackFnc);
    end

    function in_addbuttonpush(toolbarHdl, sepTag, callbackFnc, imgFil, tooltipTxt)
        if ischar(callbackFnc) || isstring(callbackFnc)
            callbackFnc = str2func(callbackFnc);
        end
        if toolbarHdl == 0
            barhandle = DeftToolbarHandle;
        elseif toolbarHdl == 1
            barhandle = MainToolbarHandle;
        elseif toolbarHdl == 2
            barhandle = UserToolbarHandle;
        end
        pt = uipushtool(barhandle, 'Separator', sepTag);
        %pt.Icon = fullfile(mfolder,'..','resources',imgFil);
        pt.CData = in_getPtImage(imgFil);
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
        pushbuttonV=[pushbuttonV; pt];
    end

    function in_addbuttontoggle(toolbarHdl, sepTag, callbackFnc)
        imgFil = callbackFnc{3};
        tooltipTxt = callbackFnc{6};
        %if ischar(callbackFnc{1}) || isstring(callbackFnc{1})
        %    callbackFnc=str2func(callbackFnc{1});
        %end
        if toolbarHdl == 0
            barhandle = DeftToolbarHandle;
        elseif toolbarHdl == 1
            barhandle = MainToolbarHandle;
        elseif toolbarHdl == 2
            barhandle = UserToolbarHandle;
        end
        pt = uitoggletool(barhandle, 'Separator', sepTag);
        pt.CData = in_getPtImage(imgFil);
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
        pushbuttonV=[pushbuttonV; pt];        
    end

    function in_togglebtfun(src, ~, func, ~, imgFil, ...
            actiondelay, tooltipTxt)
        if nargin < 6, actiondelay = true; end
        src.CData = in_getPtImage(imgFil);
        if actiondelay
            if src.State == "off"
                func(src);
            else
                s = 'To execute the function, click the button again or locate and click the same button in the toolbar above. Hover over the button to view a description of its function.';
                    uiwait(helpdlg(sprintf('%s\n%s', upper(tooltipTxt), s), ''));
            end
        else
            func(src);
        end
    end

    function [ptImage] = in_getPtImage(imgFil)
        try
            [img, map] = imread(fullfile(mfolder, 'resources', imgFil));
            ptImage = ind2rgb(img, map);
        catch
            try
                [img, map] = imread(fullfile(matlabroot,'toolbox', ...
                    'matlab','icons', imgFil));
                 ptImage = ind2rgb(img, map);
            catch
                ptImage = rand(16, 16, 3);
            end
        end
    end

    % ------------------------
    % Callback Functions
    % ------------------------

    function in_call_scgeatool(~, ~)
        % scgeatool;
        % P = get(FigureHandle,'Position');
        % k=1;
        % set(FigureHandle,'Position',[P(1)-30*k P(2)-30*k P(3) P(4)]);
    end

    function in_closeRequest(hObject, ~)
        if ~(ismcc || isdeployed)
            if isempty(sce)||sce.NumCells==0
                ButtonName='no';
            else
                ButtonName = questdlg('Save SCE before closing SCGEATOOL?');
            end
            switch lower(ButtonName)
                case 'yes'
                    if ~isempty(callinghandle)
                        guidata(callinghandle, sce);
                        delete(hObject);
                        helpdlg('SCE updated.');
                    else
                        labels = {'Save SCE to variable named:'};
                        vars = {'sce'};
                        sce = guidata(FigureHandle);
                        values = {sce};
                        [~, tf] = export2wsdlg(labels, vars, values, ...
                            'Save Data to Workspace');
                        if tf
                            delete(hObject);
                        else
                            return;
                        end
                    end
                case 'cancel'
                    return;
                case 'no'
                    delete(hObject);
                otherwise
                    return;
            end
        else
            delete(hObject);
        end
    end

    function in_SimulateSCE(src, ~)
        if ~isempty(sce) && sce.NumCells>0
            answer = questdlg('Current SCE will be replaced. Continue?');
            if ~strcmp(answer, 'Yes'), return; end
        end
        definput = {'3000', '5000'};
        prompt = {'Number of genes:', ...
            'Number of cells:'};
        dlgtitle = 'Simulation Settings';
        dims = [1, 55];
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer), return; end
        try
            numgenes = str2double(answer{1});
            numcells = str2double(answer{2});
            assert(isfinite(numgenes) & numgenes==floor(numgenes));
            assert(isfinite(numcells) & numcells==floor(numcells));
            assert((numgenes >= 1) && (numgenes <= 30000));
            assert((numcells >= 1) && (numcells <= 30000));
        catch
            errordlg('Invalid parameter value(s).');
            return;
        end
        
            try
                fw = gui.gui_waitbar;
                [X] = sc_simudata(numgenes, numcells,'lun');
                [sce] = SingleCellExperiment(X);
                [c, cL] = grp2idx(sce.c);
                gui.gui_waitbar(fw);
                guidata(FigureHandle, sce);
                in_RefreshAll(src, [], false, false);
            catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
            end

    end

    function in_GEOAccessionToSCE(src, ~)
        answer = questdlg('Current SCE will be replaced. Continue?');
        if ~strcmp(answer, 'Yes'), return; end       
        acc = inputdlg({'Input Number(s) (e.g., GSM3308547,GSM3308548):'}, ...
                    'GEO Accession', [1, 50], {'GSM3308547'});        
        if isempty(acc), return; end
        acc = acc{1};
        if strlength(acc) > 4 && ~isempty(regexp(acc, 'G.+', 'once'))
            try
                fw = gui.gui_waitbar;
                [sce] = sc_readgeoaccession(acc);
                [c, cL] = grp2idx(sce.c);
                gui.gui_waitbar(fw);
                guidata(FigureHandle, sce);
                in_RefreshAll(src, [], false, false);
            catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
            end
        end
    end

    function in_RunDataMapPlot(src, ~)        
        ndim = 2;
        [vslist] = gui.i_checkexistingembed(sce, ndim);
        if isempty(h.ZData) && size(sce.s,2)==2 && length(vslist) <= 1
            gui.callback_RunDataMapPlot(src, []);
        elseif isempty(h.ZData) && size(sce.s,2)==2 && length(vslist) > 1
            answer = questdlg('Using current 2D embedding?');
            switch answer
                case 'Yes'
                    gui.callback_RunDataMapPlot(src, []);
                case 'No'
                    [sx] = gui.i_pickembedvalues(sce, 2);
                    if ~isempty(sx) && size(sx,1) == sce.NumCells
                        sce.s = sx;
                    else
                        warning('Running error.');
                        return;
                    end
                    guidata(FigureHandle, sce);
                    gui.callback_RunDataMapPlot(src, []);
                case 'Cancel'
                    return;
            end
        elseif ~isempty(h.ZData)
            answer=questdlg('This function requires 2D embedding. Continue?');
            switch answer
                case 'Yes'
                    in_Switch2D3D(src,[]);
                otherwise
                    return;
            end
        end
    end

    function in_SubtypeAnnotation(src, ~)
        [requirerefresh] = gui.callback_SubtypeAnnotation(src, []);
        if requirerefresh
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c_cell_type_tx);
            in_RefreshAll(src, [], true, false);
            ix_labelclusters(true);
        end        
    end

    function in_MergeCellSubtypes(src, ~, sourcetag, allcell)
        if nargin < 4
            answer = questdlg('Import annotation for all cells or just cells of a subtype?', '', ...
                    'All Cells', 'Subtype Cells', 'Cancel', 'All Cells');
            switch answer
                case 'All Cells'
                    allcell = true;
                case 'Subtype Cells'
                    allcell = false;
                case 'Cancel'
                    return;
            end
        end
        [requirerefresh] = gui.callback_MergeCellSubtypes(src, [], sourcetag, allcell);
        if requirerefresh
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c_cell_type_tx);
            in_RefreshAll(src, [], true, false);
            ix_labelclusters(true);
        end
    end

    function in_MergeSCEs(src, ~, sourcetag)
        [requirerefresh, s] = gui.callback_MergeSCEs(src, sourcetag);
        if requirerefresh && ~isempty(s)
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c_batch_id);
            sce.c = c;
            in_RefreshAll(src, [], true, false);
            helpdlg(sprintf('%s SCEs merged.', upper(s)), '');
        end
    end

    function in_WorkonSelectedGenes(src, ~)
            %         answer=questdlg('Input the number of HVGs. Continue?');
            %         if ~strcmp(answer,'Yes'), return; end
            k = gui.i_inputnumk(2000, 1, sce.NumGenes, 'the number of HVGs');
            if isempty(k), return; end
            answer = questdlg('Which method?', 'Select Method', ...
                'Brennecke et al. (2013)', 'Splinefit Method', ...
                'Brennecke et al. (2013)');
            fw = gui.gui_waitbar;
            switch answer
                case 'Brennecke et al. (2013)'
                    T = sc_hvg(sce.X, sce.g);
                case 'Splinefit Method'
                    T = sc_splinefit(sce.X, sce.g);
                otherwise
                    return;
            end
            glist = T.genes(1:min([k, sce.NumGenes]));
            [y, idx] = ismember(glist, sce.g);
            if ~all(y), errordlg('Runtime error.');
                return;
            end
            sce.g = sce.g(idx);
            sce.X = sce.X(idx, :);
            gui.gui_waitbar(fw);
            in_RefreshAll(src, [], true, false);
    end

    function in_SubsampleCells(src, ~, methodoption)
        if nargin < 3
            methodoption = [];
        end
        answer1 = questdlg('This function subsamples 50% of cells. Continue?');
        if ~strcmp(answer1, 'Yes')
            return;
        end

        if isempty(methodoption)
            answer = questdlg('Select method:', '', ...
                'Uniform Sampling', ...
                'Geometric Sketching [PMID:31176620]', 'Uniform Sampling');
            switch answer
                case 'Uniform Sampling'
                    methodoption = 1;
                case 'Geometric Sketching [PMID:31176620]'
                    methodoption = 2;
                otherwise
                    return;
            end
        end

        tn = round(sce.NumCells/2);
        if methodoption == 1
            idx = randperm(sce.NumCells);
            ids = idx(1:tn);
        elseif methodoption == 2
            gui.gui_showrefinfo('Geometric Sketching [PMID:31176620]');
            fw = gui.gui_waitbar;
            Xn = log(1+sc_norm(sce.X))';
            [~, Xn] = pca(Xn, 'NumComponents', 300);
            gui.gui_waitbar(fw);
            try
                ids = run.py_geosketch(Xn, tn);
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end

        end
        if ~isempty(ids)
            sce = sce.selectcells(ids);
            c = sce.c;
            in_RefreshAll(src, [], true, false);
        else
            errordlg('Running error. No action is taken.');
        end
    end

    function in_SingleClickSolution(src, ~)
        speciestag = gui.i_selectspecies(2);
        if isempty(speciestag), return; end

        fw = gui.gui_waitbar_adv;
        gui.gui_waitbar_adv(fw,1/8,'Basic QC Filtering...');
        sce = sce.qcfilter;
        gui.gui_waitbar_adv(fw,2/8, 'Embeding Cells Using tSNE...');
        sce = sce.embedcells('tsne3d', true, true, 3);

        gui.gui_waitbar_adv(fw,3/8, 'Clustering Cells Using K-means...');
        sce = sce.clustercells([], [], true);
        gui.gui_waitbar_adv(fw,4/8, 'Annotating Cell Type Using PanglaoDB...');
        
        tic
        sce = sce.assigncelltype(speciestag, false);
        toc
        gui.gui_waitbar_adv(fw,5/8, 'Estimate Cell Cycles...');
        
        sce = sce.estimatecellcycle;
        gui.gui_waitbar_adv(fw,6/8, 'Estimate Differentiation Potency of Cells...');

        sce = sce.estimatepotency(speciestag);

        gui.gui_waitbar_adv(fw,7/8);
        [c,cL] = grp2idx(sce.c_cell_type_tx);
        sce.c = c;
        gui.gui_waitbar_adv(fw);
        in_RefreshAll(src, [], true, false);
        ix_labelclusters(true);
        %setappdata(FigureHandle, 'cL', cL);
        guidata(FigureHandle, sce);
    end

    function in_SelectCellsByQC(src, ~)
        oldn = sce.NumCells;
        oldm = sce.NumGenes;
        sce.c = c;
        guidata(FigureHandle, sce);
        try
        [requirerefresh, highlightindex] = ...
            gui.callback_SelectCellsByQC(src);
        catch ME
            errordlg(ME.message);
            return;
        end
        sce = guidata(FigureHandle);
        if requirerefresh
            [c, cL] = grp2idx(sce.c);
            in_RefreshAll(src, [], true, false);
            newn = sce.NumCells;
            newm = sce.NumGenes;
            helpdlg(sprintf('%d cells removed; %d genes removed.', ...
                oldn-newn, oldm-newm), '');
        end
        if ~isempty(highlightindex)
            h.BrushData = highlightindex;
        end
    end

    function in_Select5000Genes(src, ~)
        oldm = sce.NumGenes;
        oldn = sce.NumCells;
        [requirerefresh] = gui.callback_Select5000Genes(src);
        
        if requirerefresh
            sce = guidata(FigureHandle);
            try
                % case 'Relaxed (keep more cells/genes)'
                %     definput = {'500', '0.20', '10', '200'};
                % case 'Strigent (remove more cells/genes)'
                %     definput = {'1000', '0.15', '15', '500'};
                sce = sce.qcfilterwhitelist(1000, 0.15, 15, 500, []);
                disp('Strigent QC applied.');
                c=sce.c;
            catch ME
                warning(ME.message);
            end

            in_RefreshAll(src, [], true, false);
            newm = sce.NumGenes;
            newn = sce.NumCells;
            helpdlg(sprintf('%d cells removed; %d genes removed.', ...
                oldn-newn, oldm-newm), '');
            
            % helpdlg(sprintf('%d genes removed.', oldm-newm), '');
        end        
    end

    function in_RunSeuratWorkflow(src, ~)
        extprogname = 'R_Seurat';
        preftagname = 'externalwrkpath';
        [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);

        [ok] = gui.i_confirmscript('Run Seurat/R Workflow (Seurat)?', ...
            'R_Seurat', 'r');
        if ~ok, return; end

        [ndim] = gui.i_choose2d3d;
        if isempty(ndim), return; end
        fw = gui.gui_waitbar;
        try
            [sce] = run.r_seurat(sce, ndim, wkdir);
            [c, cL] = grp2idx(sce.c);
        catch
            gui.gui_waitbar(fw);
            return;
        end
        gui.gui_waitbar(fw);
        guidata(FigureHandle, sce);
        in_RefreshAll(src, [], true, false);
    end

    function in_DecontX(~, ~)
        gui.gui_showrefinfo('DecontX [PMID:32138770]');

        extprogname = 'R_decontX';
        preftagname = 'externalwrkpath';
        [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);

        [ok] = gui.i_confirmscript('Detect Ambient RNA Contamination (decontX)', ...
            'R_decontX', 'r');
        if ~ok, return; end
        fw = gui.gui_waitbar;
        try
            [Xdecon, contamination] = run.r_decontX(sce, wkdir);
        catch
            gui.gui_waitbar(fw);
            errordlg('Runtime error.')
            return;
        end
        gui.gui_waitbar(fw);
        figure('WindowStyle', 'modal');
        gui.i_stemscatter(sce.s, contamination);
        % zlim([0 1]);
        zlabel('Contamination rate')
        title('Ambient RNA contamination')
        answer = questdlg("Remove contamination?");
        switch answer
            case 'Yes'
                sce.X = round(Xdecon);
                guidata(FigureHandle, sce);
                helpdlg('Contamination removed.', '');
        end
    end

    function in_HarmonyPy(src, ~)
        gui.gui_showrefinfo('Harmony [PMID:31740819]');
        if numel(unique(sce.c_batch_id)) < 2
            warndlg('No batch effect (SCE.C_BATCH_ID is empty)');
            return;
        end
        [c1] = grp2idx(sce.c);
        [c2] = grp2idx(sce.c_batch_id);
        if ~isequal(c1, c2)
            answer = questdlg('Color cells by batch id (SCE.C_BATCH_ID)?', '');
            switch answer
                case 'Yes'
                    [c, cL] = grp2idx(sce.c_batch_id);
                    sce.c = c;
                    in_RefreshAll(src, [], true, false);
                case 'No'
                case 'Cancel'
                    return;
                otherwise
                    return;
            end
        end
        
        if gui.callback_Harmonypy(src)
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c);
            in_RefreshAll(src, [], true, false);
    
            ButtonName = questdlg('Update Saved Embedding?', '');
            switch ButtonName
                case 'Yes'
                    [methodtag] = gui.i_pickembedmethod;
                    if isempty(methodtag), return; end
                    [ndim] = gui.i_choose2d3d;
                    if isempty(ndim), return; end
                    methoddimtag = sprintf('%s%dd',methodtag, ndim);
                    if ismember(methoddimtag, fieldnames(sce.struct_cell_embeddings))
                        sce.struct_cell_embeddings.(methodtag) = sce.s;
                    end
                    helpdlg(sprintf('%s Embedding is updated.', methoddimtag), '');
            end
        end
        guidata(FigureHandle, sce);
    end

    function in_DoubletDetection(src, ~)
        gui.gui_showrefinfo('Scrublet [PMID:30954476]');
        [isDoublet, doubletscore, methodtag, done] = gui.callback_DoubletDetection(src);
        if done && ~any(isDoublet)
            helpdlg('No doublet detected.', '');
            return;
        end
        if done && any(isDoublet) && sce.NumCells == length(doubletscore)
            tmpf_doubletdetection = figure('WindowStyle', 'modal');
            gui.i_stemscatter(sce.s, doubletscore);
            zlabel('Doublet Score')
            title(sprintf('Doublet Detection (%s)', methodtag))
            answer = questdlg(sprintf("Remove %d doublets?", sum(isDoublet)));
            switch answer
                case 'Yes'
                    close(tmpf_doubletdetection);
                    % i_deletecells(isDoublet);
                    sce = sce.removecells(isDoublet);
                    guidata(FigureHandle, sce);
                    [c, cL] = grp2idx(sce.c);
                    in_RefreshAll(src, [], true, false);
                    helpdlg('Doublets deleted.', '');
            end
        end
    end

    function in_MergeSubCellTypes(src, ~)
        if isempty(sce.c_cell_type_tx), return; end
        newtx = erase(sce.c_cell_type_tx, "_{"+digitsPattern+"}");
        if isequal(sce.c_cell_type_tx, newtx)
            helpdlg("No sub-clusters are meraged.");
        else
            sce.c_cell_type_tx = newtx;
            [c, cL] = grp2idx(sce.c_cell_type_tx);
            sce.c = c;
            in_RefreshAll(src, [], true, false);
            ix_labelclusters(true);
        end
        guidata(FigureHandle, sce);
    end

    function in_EnDisableMenu(entag)
        % if nargin<1, entag='off'; end
        % for k=1:length(pushbuttonV)
        %     set(pushbuttonV(k),'Enable',entag);
        % end
        set(DeftToolbarHandle,'Visible',entag);
        set(MainToolbarHandle,'Visible',entag);
        showuseronboarding = getpref('scgeatoolbox', 'useronboardingtoolbar');
        switch entag
            case 'on'
                if showuseronboarding
                    set(UserToolbarHandle,'Visible','on');
                else
                    set(UserToolbarHandle,'Visible','off');
                end
            case 'off'
                set(UserToolbarHandle,'Visible','off');
        end
        menusv={m_file,m_edit,m_view,m_plot,m_ext,m_exp,m_net};
        for j=1:length(menusv)
            a=allchild(menusv{j});
            for k=1:length(a)
                a(k).Enable=entag;
            end
        end
        a=allchild(m_file);
        a(end).Enable='on';
        a(end-1).Enable='on';
        %a(end-3).Enable='on';
        a(end-6).Enable='on';
        %allchild(m_file)
        a=allchild(m_ext);
        a(end).Enable='on';
        a(end-1).Enable='on';
        a(end-2).Enable='on';

    end

    function in_RefreshAll(src, ~, keepview, keepcolr)
        if nargin < 4, keepcolr = false; end
        if nargin < 3, keepview = false; end
        if keepview || keepcolr
            [para] = gui.i_getoldsettings(src);
        end
        if isempty(sce) || sce.NumCells==0
            set(hAx,'Visible','off');
            in_EnDisableMenu('off');
            return;
        end
        in_EnDisableMenu('on');        
        figure(FigureHandle);
        % [c,cL]=grp2idx(sce.c);
        % was3d = ~isempty(h.ZData);
        if size(sce.s, 2) >= 3
            if keepview, [ax, bx] = view(hAx); end
            if isempty(c), [c,cL] = grp2idx(sce.c); end
            h = gui.i_gscatter3(sce.s, c, methodid, hAx);
            if keepview, view(ax, bx); end
        else        % otherwise going to show 2D            
            if keepview, [ax, bx] = view(hAx); end
            h = gui.i_gscatter3(sce.s(:, 1:2), c, methodid, hAx);
            if keepview, [ax, bx] = view(hAx); end
        end
        if keepview
            h.Marker = para.oldMarker;
            h.SizeData = para.oldSizeData;
        end
        if keepcolr
            colormap(para.oldColorMap);
        else
            kc = numel(unique(sce.c));
            colormap(pkg.i_mycolorlines(kc));
        end
        title(hAx, sce.title);
        subtitle(hAx, '[genes x cells]');
        set(hAx,'Visible','on');
    end

    function in_Switch2D3D(src, ~)  
        [para] = gui.i_getoldsettings(src);

        if isempty(h.ZData)               % current 2D
            ansx = questdlg('Switch to 3D?');
            if ~strcmp(ansx, 'Yes'), return; end
            if size(sce.s, 2) >= 3
                h = gui.i_gscatter3(sce.s, c, methodid, hAx);
                if ~isempty(ax) && ~isempty(bx) && ~any([ax, bx] == 0)
                    view(ax, bx);
                else
                    view(3);
                end
            else
                [vslist] = gui.i_checkexistingembed(sce, 3);
                if isempty(vslist)
                    in_EmbeddingAgain(src, [], 3);
                else
                    ansx = questdlg('Using existing 3D embedding? Select "No" to re-embed.');
                    switch ansx
                        case 'Yes'
                            [sx] = gui.i_pickembedvalues(sce, 3);
                            if ~isempty(sx) && size(sx,1) == sce.NumCells
                                sce.s = sx;
                            else
                                warning('Running error.');
                                return;
                            end
                        case 'No'
                            in_EmbeddingAgain(src, [], 3);
                        case 'Cancel'
                            return;
                    end
                end
            end
        else        % current 3D do following
            ansx = questdlg('Switch to 2D?');
            if ~strcmp(ansx, 'Yes'), return; end
            [vslist] = gui.i_checkexistingembed(sce, 2);
            if ~isempty(vslist)
                answer = questdlg('How to make 2D embedding?','', ...
                    'Pick existing 2D','Re-embed cells', ...
                    'Reduce current 3D','Pick existing 2D');
            else
                answer = questdlg('How to make 2D embedding?','', ...
                    'Re-embed cells', ...
                    'Reduce current 3D to 2D','Cancel','Re-embed cells');
            end
            switch answer
                case 'Cancel'
                    return;
                case 'Re-embed cells'
                    in_EmbeddingAgain(src, [], 2);
                case 'Reduce current 3D to 2D'
                    [ax, bx] = view(hAx);
                    answer2 = questdlg('Which view to be used to project cells?', '', ...
                        'X-Y Plane', 'Screen/Camera', 'PCA-rotated', 'X-Y Plane');
                    switch answer2
                        case 'X-Y Plane'
                            sx = sce.s;
                        case 'Screen/Camera'
                            sx = pkg.i_3d2d(sce.s, ax, bx);
                        case {'PCA-rotated'}
                            [~, sx] = pca(sce.s);
                        otherwise
                            return;
                    end
                    h = gui.i_gscatter3(sx(:, 1:2), c, methodid, hAx);
                    sce.s = sx(:, 1:2);
                    title(hAx, sce.title);
                    subtitle(hAx, '[genes x cells]');
                    h.Marker = para.oldMarker;
                    h.SizeData = para.oldSizeData;
                    colormap(para.oldColorMap);
                    return;                    
                case 'Pick existing 2D'
                    [sx] = gui.i_pickembedvalues(sce, 2);
                    if ~isempty(sx) && size(sx,1) == sce.NumCells
                        sce.s = sx;
                    else
                        warning('Running error.');
                        return;
                    end
            end
        end        
        guidata(FigureHandle, sce);
        in_RefreshAll(src, [], true, true);   % keepview, keepcolr       
    end

    function in_AddEditCellAttribs(~,~)
        answer = questdlg('Add or edit cell attribute?','', ...
            'Add','Edit','Cancel','Add');
        switch answer
            case 'Edit'
                addnew = false;
            case 'Add'
                addnew = true;
            otherwise
                return;
        end        
        [sce] = gui.sc_cellattribeditor(sce, addnew);
        guidata(FigureHandle, sce);
    end

    function in_ExportCellAttribTable(~,~)
        %sce = guidata(FigureHandle);
        T = pkg.makeattributestable(sce);        
        gui.i_exporttable(T,true, ...
            "Tcellattrib","CellAttribTable");
    end

    function in_RenameCellTypeBatchID(src, ~, answer)
        if nargin < 3 || isempty(answer)
            answer = questdlg('Rename cell type, batch ID, or gene name?', '', 'Cell type', 'Batch ID', 'Gene name', 'Cell type');
        end
        switch answer
            case 'Cell type'
                [requirerefresh] = gui.callback_RenameCellType(src);
            case 'Batch ID'
                [requirerefresh] = gui.callback_RenameBatchID(src);
            case 'Gene name'
                [requirerefresh] = gui.callback_RenameGenes(src);
            otherwise
                return;
        end
        if requirerefresh
            sce = guidata(FigureHandle);
            switch answer
                case 'Cell type'
                    [c, cL] = grp2idx(sce.c_cell_type_tx);
                case 'Batch ID'
                    [c, cL] = grp2idx(sce.c_batch_id);
                otherwise
                    return;
            end
            ix_labelclusters(false);
        end
    end

    function in_EmbeddingAgain(src, ~, ndim)
        if nargin<3, ndim=[]; end
        if isempty(sce.struct_cell_embeddings)
            sce.struct_cell_embeddings = pkg.e_makeembedstruct;
        end
        if ~isfield(sce.struct_cell_embeddings, 'metaviz3d')
            sce.struct_cell_embeddings.('metaviz3d') = [];
        end

        [vslist] = gui.i_checkexistingembed(sce, ndim);
        if ~isempty(vslist)
            answer = questdlg('Using exsiting embedding?');
        else
            answer = 'No';
        end
        if strcmp(answer, 'Yes')
            [sx] = gui.i_pickembedvalues(sce);
            if ~isempty(sx), sce.s = sx; end
        elseif strcmp(answer, 'No')
            [methodtag] = gui.i_pickembedmethod;
            if isempty(methodtag), return; end
            if isempty(ndim), [ndim] = gui.i_choose2d3dnmore; end
            if isempty(ndim), return; end
            methoddimtag = sprintf('%s%dd',methodtag, ndim);
            usingold = false;
            if ~isfield(sce.struct_cell_embeddings, methoddimtag)
                sce.struct_cell_embeddings = setfield(sce.struct_cell_embeddings,methoddimtag,[]);
            end
            
            % if ~isempty(sce.struct_cell_embeddings.(methoddimtag))
            %     answer1 = questdlg(sprintf('Use existing %s embedding or re-compute new embedding?', ...
            %         upper(methoddimtag)), '', ...
            %         'Use existing', 'Re-compute', 'Cancel', 'Use existing');
            %     switch answer1
            %         case 'Use existing'
            %             sce.s = sce.struct_cell_embeddings.(methoddimtag);
            %             usingold = true;
            %         case 'Re-compute'
            %             usingold = false;
            %         case {'Cancel', ''}
            %             return;
            %     end
            % end

            % whitelist = [];
            if ~usingold
                [K, usehvgs] = gui.i_gethvgnum(sce);
                if isempty(K), return; end
                fw = gui.gui_waitbar;
                try
                    forced = true;
                    %if contains(methoddimtag, 'tsne'), disp('tSNE perplexity = 30'); end
                    sce = sce.embedcells(methodtag, forced, usehvgs, ndim, K);
                    % disp('Following the library-size normalization and log1p-transformation, we visualized similarity among cells by projecting them into a reduced dimensional space using t-distributed stochastic neighbor embedding (t-SNE)/uniform manifold approximation and projection (UMAP).')
                catch ME
                    gui.gui_waitbar(fw, true);
                    errordlg(ME.message);
                    % rethrow(ME)
                    return;
                end
                gui.gui_waitbar(fw);
            end
        else
            return;
        end
        guidata(FigureHandle, sce);
        in_RefreshAll(src, [], true, false);   % keepview, keepcolr
    end

    function in_DetermineCellTypeClustersGeneral(src, ~, usedefaultdb)
        if nargin < 3, usedefaultdb = true; end
        if usedefaultdb
            organtag = "all";
            databasetag = "panglaodb";
            gui.gui_showrefinfo('PanglaoDB [PMID:30951143]');
            speciestag = gui.i_selectspecies(2);
            if isempty(speciestag), return; end
        else
            [Tm, Tw] = pkg.i_markerlist2weight(sce);
            if isempty(Tm) || isempty(Tw)
                return;
            end
            wvalu = Tw.Var2;
            wgene = string(Tw.Var1);
            celltypev = string(Tm.Var1);
            markergenev = string(Tm.Var2);
        end

        [manuallyselect, bestonly] = i_annotemanner;

        dtp = findobj(h, 'Type', 'datatip');
        delete(dtp);
        cLdisp = cL;
        if ~manuallyselect, fw = gui.gui_waitbar_adv; end

        for ix = 1:max(c)
            if ~manuallyselect
                gui.gui_waitbar_adv(fw, ix/max(c));
            end
            ptsSelected = c == ix;

            if usedefaultdb
                [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, ...
                    sce.s, ptsSelected, ...
                    speciestag, organtag, databasetag, bestonly);
            else
                [Tct] = pkg.e_determinecelltype(sce, ptsSelected, wvalu, ...
                    wgene, celltypev, markergenev);
            end

            ctxt = Tct.C1_Cell_Type;

            if manuallyselect && length(ctxt) > 1
                [indx, tf] = listdlg('PromptString', {'Select cell type'}, ...
                    'SelectionMode', 'single', 'ListString', ctxt);
                if tf ~= 1, return; end
                ctxt = Tct.C1_Cell_Type{indx};
            else
                ctxt = Tct.C1_Cell_Type{1};
            end

            hold on;
            ctxtdisp = strrep(ctxt, '_', '\_');
            ctxtdisp = sprintf('%s_{%d}', ctxtdisp, ix);
            cLdisp{ix} = ctxtdisp;

            ctxt = sprintf('%s_{%d}', ctxt, ix);
            cL{ix} = ctxt;

            row = dataTipTextRow('', cLdisp(c));
            h.DataTipTemplate.DataTipRows = row;
            if size(sce.s, 2) >= 2
                siv = sce.s(ptsSelected, :);
                si = mean(siv, 1);
                idx = find(ptsSelected);
                [k] = dsearchn(siv, si); % Nearest point search
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
        if ~manuallyselect, gui.gui_waitbar_adv(fw); end
        sce.c_cell_type_tx = string(cL(c));
        nx = length(unique(sce.c_cell_type_tx));
        if nx > 1
            newtx = erase(sce.c_cell_type_tx, "_{"+digitsPattern+"}");
            if length(unique(newtx)) ~= nx
                answer = questdlg('Merge subclusters of same cell type?');
                if strcmp(answer, 'Yes')
                    in_MergeSubCellTypes(src);
                end
            end
        end
        guidata(FigureHandle, sce);
    end

    function [iscelltype] = in_pickcelltypeclusterid(a)
        answer = questdlg(a, '', ...
            'Cluster', 'Cell Type', 'Cluster');
        switch answer
            case 'Cluster'
                iscelltype = false;
            case 'Cell Type'
                iscelltype = true;
            otherwise
                iscelltype = [];
        end
    end

    function in_Brushed2NewCluster(~, ~)
        [iscelltype] = in_pickcelltypeclusterid('Make a new cluster or new cell type group out of brushed cells?');
        if isempty(iscelltype), return; end
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return;
        end
        if iscelltype
            n = sum(contains(unique(sce.c_cell_type_tx), "New cell type"));
            if n > 0
                nname = sprintf('New cell type %d', n+1);
            else
                nname = 'New cell type';
            end
            newctype = inputdlg('Enter new cell type name:', '', [1, 50], {nname});
            if isempty(newctype), return; end
            sce.c_cell_type_tx(ptsSelected) = string(newctype);
            [c, cL] = grp2idx(sce.c_cell_type_tx);
        else
            c(ptsSelected) = max(c) + 1;
            [c, cL] = grp2idx(c);
            sce.c_cluster_id = c;
        end
        sce.c = c;
        [ax, bx] = view(hAx);
        [h] = gui.i_gscatter3(sce.s, c, methodid, hAx);
        title(hAx, sce.title);
        subtitle(hAx, '[genes x cells]');
        view(ax, bx);
        ix_labelclusters(true);
        guidata(FigureHandle, sce);
    end

    function in_Brushed2MergeClusters(~, ~)
        [iscelltype] = in_pickcelltypeclusterid('Merge brushed cells into same cluster or same cell type?');
        if isempty(iscelltype), return; end
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are brushed");
            return;
        end
        if iscelltype
            c_members = unique(sce.c_cell_type_tx(ptsSelected));
        else
            c_members = unique(c(ptsSelected));
        end

        if numel(c_members) == 1
            warndlg("All brushed cells are in one cluster or belong to the same cell type.");
            return;
        end

        [indx, tf] = listdlg('PromptString', ...
            {'Select target cluster'}, 'SelectionMode', ...
            'single', 'ListString', string(c_members));
        if tf == 1
            c_target = c_members(indx);
        else
            return;
        end

        if iscelltype
            sce.c_cell_type_tx(ismember(sce.c_cell_type_tx, c_members)) = c_target;
            [c, cL] = grp2idx(sce.c_cell_type_tx);
        else
            c(ismember(c, c_members)) = c_target;
            [c, cL] = grp2idx(c);
            sce.c = c;
            sce.c_cluster_id = c;
        end
        [ax, bx] = view(hAx);
        [h] = gui.i_gscatter3(sce.s, c, methodid, hAx);
        title(hAx, sce.title);
        subtitle(hAx, '[genes x cells]');
        view(ax, bx);
        ix_labelclusters(true);
        guidata(FigureHandle, sce);
    end


    function in_Brush4Celltypes(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            helpdlg("No cells are selected. Please use the data brush tool to select cells for cell type assignment.", '');
                return;
        end
        answer = questdlg('This is a one-time analysis. Cell type labels will not be saved. Continue?');
        if ~strcmp(answer, 'Yes')
            return;
        end
        speciestag = gui.i_selectspecies;
        if isempty(speciestag), return; end
        fw = gui.gui_waitbar;
        [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, sce.s, ...
            ptsSelected, ...
            speciestag, "all", "panglaodb", false);
        ctxt = Tct.C1_Cell_Type;
        gui.gui_waitbar(fw);

        [indx, tf] = listdlg('PromptString', ...
            {'Select cell type'}, 'SelectionMode', 'single', 'ListString', ctxt);
        if tf == 1
            ctxt = Tct.C1_Cell_Type{indx};
        else
            return;
        end
        ctxt = strrep(ctxt, '_', '\_');
        delete(findall(FigureHandle, 'Type', 'hggroup'));
        if ~exist('tmpcelltypev', 'var')
            tmpcelltypev = cell(sce.NumCells, 1);
        end
        siv = sce.s(ptsSelected, :);
        si = mean(sce.s(ptsSelected, :));
        [k] = dsearchn(siv, si);
        idx = find(ptsSelected);
        tmpcelltypev{idx(k)} = ctxt;
        row = dataTipTextRow('', tmpcelltypev);
        h.DataTipTemplate.DataTipRows = row;
        datatip(h, 'DataIndex', idx(k));
    end

    function in_ShowCellStates(src, ~)
        % sce = guidata(FigureHandle);
        [thisc, clable, ~, newpickclable] = gui.i_select1state(sce);
        %clable
        %newpickclable
        if strcmp(clable, 'Cell Cycle Phase')
            if ~all(strcmp(unique(thisc), "undetermined"))
                sce.c_cell_cycle_tx = thisc;
            end
        end
        if isempty(thisc), return; end
        if strcmp(clable, 'Workspace Variable...')
            clable = gui.i_renamec(clable, sce, newpickclable);
            sce.list_cell_attributes = [sce.list_cell_attributes, clable];
            sce.list_cell_attributes = [sce.list_cell_attributes, thisc];
        end
        [c, cL] = grp2idx(thisc);
        sce.c = c;
        answer1 = questdlg('Display in place or in new figure?', '', ...
            'In place', 'New figure','Cancel','In place');
        switch answer1
            case 'In place'
                in_RefreshAll(src, [], true, false);
                target{1} = hAx;
                target{2} = h;
            case 'New figure'
                hFig = figure('Visible','off');
                hFigAx = axes('Parent', hFig);
                h2 = gui.i_gscatter3(sce.s, c, 1, 1, hFigAx);
                title(hFigAx, sce.title);
                subtitle(hFigAx, '[genes x cells]');
                target{1} = hFigAx;
                target{2} = h2;
                dp = get(hFig, 'Position');
                pp = get(FigureHandle, 'Position');
                cpx = pp(1) + pp(3)/2 - dp(3)/2;
                cpy = pp(2) + pp(4)/2 - dp(4)/2;
                movegui(hFig, [cpx cpy]);
                set(hFig,'Visible','on');
            otherwise
                return;
        end
        [answer] = gui.i_selvariabletype(thisc);        
        switch answer
            case 'Categorical/Discrete'
                n = max(c);
                f = 0.5 * (n - 1) ./ n;
                f = 1 + f .* (1:2:2 * n);
                cb = colorbar(target{1}, 'Ticks', f, 'TickLabels', ...
                    strrep(cellstr(cL), '_', '\_'));
            otherwise  % case 'Numerical/Continuous'
                if isnumeric(thisc)
                    set(target{2}, 'CData', thisc);
                else
                    set(target{2}, 'CData', c);
                end
                cb = colorbar(target{1});
        end
        cb.Label.String = strrep(clable, '_', '\_');
        guidata(FigureHandle, sce);
    end

    function in_EnrichrHVGs(src, events)
        gui.gui_showrefinfo('HVG Functional Analysis [PMID:31861624]');
        
        answer = questdlg('This function applies to a homogeneous group of cells. Remove lowly expressed genes before applying. Continue?');
        if ~strcmp(answer, 'Yes'), return; end
        
        ptsSelected = logical(h.BrushData.');
        if any(ptsSelected)
            [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, sce);
            if ~letdoit, return; end
            if sum(ptsSelected) < 200
                answer = questdlg(sprintf('Too few cells (n = %d) selected, continue?', sum(ptsSelected)));
                if ~strcmp(answer, 'Yes'), return; end
            end           
            scetmp = sce.removecells(~ptsSelected);
            scetmp = scetmp.qcfilter(1000, 0.15, 15);
            gui.callback_EnrichrHVGs(src, events, scetmp);
        else
            answer = questdlg(sprintf('All cells (n = %d) included, continue?', sce.NumCells));
            if ~strcmp(answer, 'Yes'), return; end
            gui.callback_EnrichrHVGs(src, events);
        end
    end

    function in_DeleteSelectedCells(src, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.", '');
            return;
        end
        [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, sce);
        if ~letdoit, return; end

        answer = questdlg('Delete selected or unselected cells?', '', ...
            'Selected', 'Unselected', 'Selected');
        if isempty(answer), return; end
        if strcmp(answer, 'Unselected')
            in_deletecells(src, ~ptsSelected);
        elseif strcmp(answer, 'Selected')
            in_deletecells(src, ptsSelected);
        else
            return;
        end
        guidata(FigureHandle, sce);
    end

    function in_deletecells(src, ptsSelected)
        needprogressbar = false;
        if sce.NumCells > 8000, needprogressbar = true; end
        if needprogressbar
            fw = gui.gui_waitbar;
        end
        sce = sce.removecells(ptsSelected);
        if needprogressbar
            gui.gui_waitbar(fw);
        end
        [c, cL] = grp2idx(sce.c);
        in_RefreshAll(src, [], true, true);
    end

    function in_scTenifoldNet(src,events,methodtag)
        if numel(unique(sce.c_cell_type_tx))>1
            answer=questdlg('This analysis is cell type-specific; however, current SCE contains multiple cell types. Continue?');
            if ~strcmp(answer,'Yes'), return; end
        end
        answer=questdlg('Subsample cells?','','Yes 🐢','No 🐇','Cancel','No 🐇');
        switch answer
            case 'No 🐇'
                if methodtag==1
                    gui.callback_scPCNet1(src,events);
                elseif methodtag==2
                    gui.callback_scTenifoldNet2lite(src,events);
                end
            case 'Yes 🐢'
                if methodtag==1
                    gui.callback_scTenifoldNet1(src,events);
                elseif methodtag==2
                    gui.callback_scTenifoldNet2(src,events);
                end
        end
    end

    function in_DrawKNNNetwork(~, ~)
        k = gui.i_inputnumk(3);
        if isempty(k), return; end
        fw = gui.gui_waitbar;
        set(0, 'CurrentFigure', FigureHandle);
        figure('WindowStyle', 'modal');
        sc_knngraph(sce.s, k, true);
        gui.gui_waitbar(fw);
    end

    function in_DrawTrajectory(src, ~)
        % waitfor(warndlg('This function should not be applied to tSNE and UMAP embeddings, as they "encourage a representation of the data as disjoint clusters, which is less meaningful for modeling continuous developmental trajectories" [PMID:25664528].', ''));
        if ~isempty(f_traj)
            answer = questdlg('Remove existing trajectory curve?');
            switch answer
                case 'Yes'
                    in_RefreshAll(src, [], true, true);  % keepview, keepcolr
                case 'No'
                otherwise
                    return;
            end
        end
        answer = questdlg('Which method?', '', 'splinefit', 'princurve', ...
            'Cancel', 'splinefit');
        switch answer
            case 'splinefit'
                dim = 1;
                [t, xyz1] = pkg.i_pseudotime_by_splinefit(sce.s, dim, false);
                pseudotimemethod = 'splinefit';
            case 'princurve'
                [t, xyz1] = pkg.i_pseudotime_by_princurve(sce.s, false);
                pseudotimemethod = 'princurve';
            otherwise
                return;
        end
        hold on;
        if size(xyz1, 2) >= 3
            f_traj = plot3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-r', 'linewidth', 2);
            text(xyz1(1, 1), xyz1(1, 2), xyz1(1, 3), 'Start', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            text(xyz1(end, 1), xyz1(end, 2), xyz1(end, 3), 'End', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        elseif size(xyz1, 2) == 2
            f_traj = plot(xyz1(:, 1), xyz1(:, 2), '-r', 'linewidth', 2);
            text(xyz1(1, 1), xyz1(1, 2), 'Start', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            text(xyz1(end, 1), xyz1(end, 2), 'End', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        end
        hold off;
        % pseudotimemethod
        % answer = questdlg('Save/Update pseudotime T in SCE', ...
        %     'Save Pseudotime', ...
        %     'Yes', 'No', 'Yes');

        tag = sprintf('%s_pseudotime', pseudotimemethod);
        % iscellstr(sce.list_cell_attributes(1:2:end))
        try
            idx = find(contains(sce.list_cell_attributes(1:2:end), tag));
        catch ME
            idx = [];
            warning(ME.message);
        end
        if ~isempty(idx)
            sce.list_cell_attributes{idx*2} = t;
            fprintf('%s is updated.\n', upper(tag));
        else
            sce.list_cell_attributes{end+1} = tag;
            sce.list_cell_attributes{end+1} = t;
            fprintf('%s is saved.\n', upper(tag));
        end
        guidata(FigureHandle, sce);

        answer = questdlg('View expression of selected genes', ...
            'Pseudotime Function', ...
            'Yes', 'No', 'Yes');
        switch answer
            case 'Yes'
                gui.sc_pseudotimegenes(sce, t);
            case 'No'
                return;
        end            
    end

    function in_ClusterCellsS(src, ~)
        answer = questdlg('Cluster cells using embedding S?');
        if ~strcmp(answer, 'Yes'), return; end

        [sx] = gui.i_pickembedvalues(sce);
        if isempty(sx), return; end

        answer = questdlg('Which method?', 'Select Algorithm', ...
            'K-means 🐇', 'SnnDpc [DOI:10.1016/j.ins.2018.03.031] 🐢', 'K-means 🐇');
        if strcmpi(answer, 'K-means 🐇')
            methodtag = "kmeans";
        elseif strcmpi(answer, 'SnnDpc [DOI:10.1016/j.ins.2018.03.031] 🐢')
            methodtag = "snndpc";
            gui.gui_showrefinfo('SnnDpc [DOI:10.1016/j.ins.2018.03.031]');
        else
            return;
        end
        in_reclustercells(src, methodtag, sx);        
        guidata(FigureHandle, sce);
    end

    function in_ClusterCellsX(src, ~)
        answer = questdlg('Cluster cells using expression matrix X?');
        if ~strcmp(answer, 'Yes'), return; end
        % methodtagvx = {'specter (31 secs) 🐇', 'sc3 (77 secs) 🐇', ...
        %     'simlr (400 secs) 🐢', ...
        %     'soptsc (1,182 secs) 🐢🐢', 'sinnlrr (8,307 secs) 🐢🐢🐢',};
        % methodtagv = {'specter', 'sc3', 'simlr', 'soptsc', 'sinnlrr'};
        methodtagvx = {'SC3 [PMID:28346451] 🐢🐢'};
        methodtagv = {'sc3'};
        [indx, tf] = listdlg('PromptString', ...
            {'Select a clustering algorithm'}, ...
            'SelectionMode', 'single', ...
            'ListString', methodtagvx);
        if tf == 1
            methodtag = methodtagv{indx};
        else
            return;
        end
        if (ismcc || isdeployed)
            if strcmp(methodtag, 'sc3')
                warndlg('SC3 is not working in standalone application.', '');
                return;
            end
        end
        in_reclustercells(src, methodtag, []);
        guidata(FigureHandle, sce);
    end

    function in_reclustercells(src, methodtag, sx)
        if nargin<3, sx = []; end
        methodtag = lower(methodtag);
        usingold = false;
        if ~isempty(sce.struct_cell_clusterings.(methodtag))
            answer1 = questdlg(sprintf('Using existing %s clustering?', ...
                upper(methodtag)), '', ...
                'Yes, use existing', 'No, re-compute', ...
                'Cancel', 'Yes, use existing');
            switch answer1
                case 'Yes, use existing'
                    sce.c_cluster_id = sce.struct_cell_clusterings.(methodtag);
                    usingold = true;
                case 'No, re-compute'
                    usingold = false;
                case 'Cancel'
                    return;
            end
        end
        if ~usingold
            defv = round(sce.NumCells/100, -1);
            defva = min([2, round(sce.NumCells/100, -2), round(sce.NumCells/20, -1)]);
            if defva == 0, defva = min([10, defv]); end
            defvb = max([round(sce.NumCells/20, -2), round(sce.NumCells/20, -1)]);
            k = gui.i_inputnumk(defv, defva, defvb);
            if isempty(k), return; end
            fw = gui.gui_waitbar;
            try
                % [sce.c_cluster_id]=sc_cluster_x(sce.X,k,'type',methodtag);
                sce = sce.clustercells(k, methodtag, true, sx);
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return
            end
            gui.gui_waitbar(fw);
        end
        [c, cL] = grp2idx(sce.c_cluster_id);
        sce.c = c;
        in_RefreshAll(src, [], true, false);
        guidata(FigureHandle, sce);
    end

    function in_labelcellgroups(src, ~)
        state = src.State;
        dtp = findobj(h, 'Type', 'datatip');
        %disp('...state...')
        if strcmp(state, 'off') || ~isempty(dtp) % switch from on to off
            %dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);
            set(src, 'State', 'off');
        else
            % sce = guidata(FigureHandle);
            [thisc, clable] = gui.i_select1class(sce,true);            
            if isempty(thisc)
                set(src, 'State', 'off');
                return;
            end
            [c, cL] = grp2idx(thisc);
            sce.c = c;
            in_RefreshAll(src, [], true, false);
            fprintf('Cells are colored by %s.\n', lower(clable));
            if max(c) <= 200
                if ix_labelclusters(true)
                    set(src, 'State', 'on');
                else
                    set(src, 'State', 'off');
                end
            else
                set(src, 'State', 'off');
                warndlg('Labels are not showing. Too many categories (n>200).');
            end
            %setappdata(FigureHandle, 'cL', cL);
            guidata(FigureHandle, sce);
            % colormap(lines(min([256 numel(unique(sce.c))])));
        end
    end

    function [txt] = i_myupdatefcnx(~, event_obj)
        % pos = event_obj.Position;
        idx = event_obj.DataIndex;
        txt = cL(c(idx));
    end

    function [isdone] = ix_labelclusters(notasking)
        if nargin < 1, notasking = true; end
        isdone = false;
        if ~isempty(cL)
            if notasking
                %stxtyes = c;
                stxtyes = cL(c);
            else
                [~, cLx] = grp2idx(c);
                if isequal(cL, cLx)
                    stxtyes = c;
                else
                    answer = questdlg(sprintf('Label %d groups with index or text?', ...
                        numel(cL)), 'Select Format', 'Index', ...
                        'Text', 'Cancel', 'Text');
                    switch answer
                        case 'Text'
                            stxtyes = cL(c);
                        case 'Index'
                            stxtyes = c;
                        otherwise
                            return;
                    end
                end
            end
            dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);
            if isstring(stxtyes) || iscellstr(stxtyes)
                % stxtyes = strrep(stxtyes, "_", "\_");
                stxtyes = strtrim(stxtyes);
            end
            row = dataTipTextRow('', stxtyes);
            h.DataTipTemplate.DataTipRows = row;
            % h.DataTipTemplate.FontSize = 5;
            for ik = 1:max(c)
                idx = find(c == ik);
                siv = sce.s(idx, :);
                si = median(siv, 1);                
                % si=geometric_median(siv');
                [kb] = dsearchn(siv, si);
                %[~, k] = medoid(siv);  geometric_median
                datatip(h, 'DataIndex', idx(kb));
            end
            %ptlabelclusters.State = 'on';
            isdone = true;
        end
    end

end