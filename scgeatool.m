function varargout = scgeatool(sce, varargin)

if usejava('jvm') && ~feature('ShowFigureWindows')
    error('MATLAB is in a text mode. This function requires a GUI-mode.');
end
if ~gui.i_installed('stats'), return; end

persistent speciestag

ptimgidx = 1;
ptImgCell = cell(2, 50);

import pkg.*
import gui.*

fx = [];
if nargin < 1
    try
        fxfun = @gui.sc_splashscreen;
        [fx, v1] = fxfun();
    catch
        fxfun = @gui.sc_simplesplash;
        [fx, v1] = fxfun();
    end
    sce = SingleCellExperiment;
else
    if isempty(sce), sce = SingleCellExperiment; end
    if ~isa(sce, 'SingleCellExperiment')
        error('requires >> sce = SingleCellExperiment(); scgeatool(sce)');
    end
    v1 = pkg.i_get_versionnum;
end

mfolder = fileparts(mfilename('fullpath'));

p = inputParser;
checkCS = @(x) isempty(x) | size(sce.X, 2) == length(x);
addRequired(p, 'sce', @(x) isa(x, 'SingleCellExperiment'));
addOptional(p, 'c', sce.c, checkCS);
addOptional(p, 's', [], checkCS);
addOptional(p, 'methodid', 1, @isnumeric);
addOptional(p, 'legacy', false, @islogical);
addOptional(p, 'useuifig', false, @islogical);
addOptional(p, 'callinghandle', []);
parse(p, sce, varargin{:});
callinghandle = p.Results.callinghandle;

c_in = p.Results.c;
s_in = p.Results.s;

methodid = p.Results.methodid;
legacy = p.Results.legacy;
useuifig = p.Results.useuifig; 

if legacy
    ptImgFile = fullfile(mfolder, 'resources', 'Images', 'ptImgFile_legacy.mat');
else
    ptImgFile = fullfile(mfolder, 'resources', 'Images', 'ptImgFile.mat');
end
if exist(ptImgFile, 'file'), load(ptImgFile, 'ptImgCell'); end

f_traj = [];   % trajectory curve
ax = []; bx = [];

tmpcelltypev = cell(sce.NumCells, 1);

if ~isempty(c_in), sce.c = c_in; end
if ~isempty(s_in), sce.s = s_in; end
[c, cL] = grp2idx(sce.c);


[FigureHandle, hAx] = gui.gui_createmainfigure(v1, useuifig);

if ~isempty(fx) && isvalid(fx), fxfun(fx,0.2); end
[button1, button2] = gui.gui_createbuttons(FigureHandle, @in_sc_openscedlg);

m_file = createMenus(FigureHandle, 1);
m_edit = createMenus(FigureHandle, 2);
m_view = createMenus(FigureHandle, 3);
m_plot = createMenus(FigureHandle, 4);
m_anno = createMenus(FigureHandle, 5);
m_tool = createMenus(FigureHandle, 6);
m_ntwk = createMenus(FigureHandle, 7);
m_extn = createMenus(FigureHandle, 8);
% m_optn = createMenus(FigureHandle, 9);

if ~isempty(fx) && isvalid(fx), fxfun(fx, 0.4); end

% axtoolbar(hAx, 'default');
if ~isempty(sce) && sce.NumCells>0
    h = gui.i_gscatter3(sce.s, c, methodid, [], hAx);    
    % h = scatter3(hAx,rand(1,100),rand(1,100),rand(1,100));
    title(hAx, sce.title);
    if sce.s>2
        rotate3d(hAx,'on');
    else
        brush(hAx,'on');
    end
    subtitle(hAx,'[genes x cells]');    
else
    h = [];
    hAx.Toolbar.Visible = 'off';
end

DeftToolbarHandle = uitoolbar('Parent', FigureHandle);
MainToolbarHandle = uitoolbar('Parent', FigureHandle);
UserToolbarHandle = uitoolbar('Parent', FigureHandle);
createPushButtons(FigureHandle);


if ~isempty(fx) && isvalid(fx), fxfun(fx, 0.6); end
pause(0.5);
if ~isempty(fx) && isvalid(fx), fxfun(fx, 0.8); end
if ~isempty(c)
    kc = numel(unique(c));
    colormap(hAx, pkg.i_mycolorlines(kc));
end

if ~isempty(sce) && sce.NumCells>0
    in_EnDisableMenu('on');
else
    in_EnDisableMenu('off');
end

if ~isempty(sce) && sce.NumCells>0, hAx.Visible="on"; end

createMenus(FigureHandle, 10);

if ~isempty(fx) && isvalid(fx), fxfun(fx, 1.0); end
pause(1);
if ~isempty(fx) && isvalid(fx), set(fx, 'visible', 'off'); end
pause(0.2);
delete(fx);
set(FigureHandle, 'visible', 'on');
if gui.i_isuifig(FigureHandle), focus(FigureHandle); end

if ~gui.i_isuifig(FigureHandle)
    uicontrol(button1);
end

% if useuifig
%     UI.Snackbar(FigureHandle,'Click Import Data or Press i to Start', ...
%         'type','normal','theme','success');
% end

guidata(FigureHandle, sce);
set(FigureHandle, 'CloseRequestFcn', @in_closeRequest);
if nargout > 0, varargout{1} = FigureHandle; end

if ~ispref('scgeatoolbox', 'useronboardingtoolbar')
    gui.gui_userguidingpref(true);
    setpref('scgeatoolbox', 'useronboardingtoolbar', true);
end
showuseronboarding = getpref('scgeatoolbox', 'useronboardingtoolbar', false);
if ~showuseronboarding, set(UserToolbarHandle, 'Visible', 'off'); end
if ~exist(ptImgFile, 'file'), save(ptImgFile, 'ptImgCell'); end

% ------------------------
% GUI Making Functions
% ------------------------
   
    function menus = createMenus(FigureHandle, id)
        switch id
            case 1
                menus = uimenu(FigureHandle, 'Text', '&File');
                in_addmenu(menus, 0, @in_sc_openscedlg, '&Import Data... ','I');
                in_addmenu(menus, 0, @in_closeRequest, '&Close','W');
                in_addmenu(menus, 1, {@gui.i_savemainfig, 3}, 'Save Figure to PowerPoint File...');
                in_addmenu(menus, 0, {@gui.i_savemainfig, 2}, 'Save Figure as Graphic File...');
                in_addmenu(menus, 0, {@gui.i_savemainfig, 1}, 'Save Figure as SVG File...');
                in_addmenu(menus, 1, @gui.callback_SaveX, 'Export && &Save Data...', 'S');
            case 2
                menus = uimenu(FigureHandle, 'Text', '&Edit');
                in_addmenu(menus, 0, @in_SelectCellsByQC, 'Filter Genes & Cells...', 'F');
                in_addmenu(menus, 1, @in_Brushed2NewCluster, 'Add Brushed Cells to a New Group');
                in_addmenu(menus, 0, @in_Brushed2MergeClusters, 'Merge Brushed Cells to Same Group');
                in_addmenu(menus, 0, @in_RenameCellTypeBatchID, 'Rename Cell Type or Batch ID...');
                in_addmenu(menus, 1, @gui.callback_SplitAtacGex, 'Split Multiome ATAC+GEX Matrix...');
                in_addmenu(menus, 1, {@in_MergeSCEs, 1}, 'Merge SCE Data Variables in Workspace...');
                in_addmenu(menus, 0, {@in_MergeSCEs, 2}, 'Merge SCE Data Files...');
                in_addmenu(menus, 1, @in_AddEditCellAttribs, 'Add/Edit Cell Attributes...');
                in_addmenu(menus, 0, @in_ExportCellAttribTable, 'Export Cell Attribute Table...');
                in_addmenu(menus, 1, @gui.callback_SelectCellsByMarker, 'Extract Cells by Marker (+/-) Expression...');
                in_addmenu(menus, 0, @in_MergeSubCellTypes, 'Merge Subclusters of Same Cell Type');
                in_addmenu(menus, 1, {@in_WorkonSelectedGenes, 'name'}, 'Select Genes to Work on...');
                in_addmenu(menus, 0, {@in_WorkonSelectedGenes, 'hvg'}, 'Select Highly Variable Genes (HVGs) to Work on...');
                in_addmenu(menus, 0, {@in_WorkonSelectedGenes, 'ligandreceptor'}, 'Select Ligand Receptor Genes to Work on...');
                in_addmenu(menus, 0, @in_SubsampleCells, 'Subsample 50% Cells to Work on...');
                in_addmenu(menus, 1, {@in_DeleteBrushedOrUnbrushedCells, 'brushed'}, 'Delete Brushed Cells...');
                in_addmenu(menus, 0, {@in_DeleteBrushedOrUnbrushedCells, 'unbrushed'}, 'Delete Unbrushed Cells...');
                in_addmenu(menus, 1, @gui.callback_SelectCellsByClass, 'Select Cells by Class & Open in New Window...');
            case 3
                menus = uimenu(FigureHandle, 'Text', '&View');
                in_addmenu(menus, 0, @in_EmbeddingAgain, 'Embed Cells Using tSNE, UMP, PHATE...', 'B');
                in_addmenu(menus, 0, @in_Switch2D3D, 'Switch Between 2D/3D Embeddings...');
                in_addmenu(menus, 1, @in_ClusterCellsS, "Cluster Cells Using Cell Embedding (S)", 'C');
                in_addmenu(menus, 0, @in_ClusterCellsX, "Cluster Cells Using Expression Matrix (X) 🐢 ...");
                in_addmenu(menus, 1, @gui.callback_ShowGeneExpr, 'Gene Expression...', 'E');
                in_addmenu(menus, 0, @in_ShowCellStates, 'Show Cell States...', 'T');
                in_addmenu(menus, 0, @in_labelcellgroups, 'Label Cell Groups...', 'G');
                in_addmenu(menus, 0, @in_highlightcellgroups, 'Highlight Cell Groups...');
                in_addmenu(menus, 0, @gui.callback_MultiGroupingView, 'Multi-Grouping View...');
                in_addmenu(menus, 0, @gui.callback_CrossTabulation, 'Cross Tabulation');
                in_addmenu(menus, 0, @gui.callback_ShowGeneExprGroup, 'Gene Expression in Groups');
                in_addmenu(menus, 1, @gui.callback_ViewMetaData, 'View Metadata...', 'M');
                in_addmenu(menus, 1, @gui.callback_ShowHgBGeneExpression, 'Hemoglobin (Hgb) Genes Expression...');
                in_addmenu(menus, 0, @gui.callback_ShowMtGeneExpression, 'Mitochondrial (Mt-) Genes Expression...');
                in_addmenu(menus, 0, @in_qcviolin, 'Cell QC Metrics in Violin Plots...');
                in_addmenu(menus, 1, @gui.callback_ShowClustersPop,"Show Cell Clusters/Groups Individually...");
                in_addmenu(menus, 1, @gui.callback_CloseAllOthers, 'Close All Other Figures', 'X');
                in_addmenu(menus, 0, @in_RefreshAll, 'Refresh Current View', 'R');
            case 4
                menus = uimenu(FigureHandle, 'Text', '&Plots');
                in_addmenu(menus, 0, @gui.callback_Dotplot, 'Gene Expression Dot Plot...');
                in_addmenu(menus, 0, @gui.callback_Heatmap, 'Gene Expression Heatmap...');
                in_addmenu(menus, 0, @gui.callback_ScatterStemPlot,'Gene Expression/Cell State Stem Plot...');
                in_addmenu(menus, 0, @gui.callback_Violinplot, 'Gene Expression/Cell State Violin Plot...');
                in_addmenu(menus, 0, @gui.callback_ScatterCorrPlot,'Correlation Plot...');
                in_addmenu(menus, 1, @gui.callback_ShowGeneExprCompr,'Side-by-Side Gene Expression...');
                in_addmenu(menus, 0, @gui.callback_EnrichrTab2Circos,'Enrichr Result Table to Circos Plot...');
                in_addmenu(menus, 1, @gui.callback_GetCellSignatureMatrix, 'Cell State Radar Plot...');
                in_addmenu(menus, 0, @in_DrawKNNNetwork, 'Cell kNN Network...');
                in_addmenu(menus, 0, @in_DrawTrajectory, 'Cell Trajectory...');
                in_addmenu(menus, 1, @gui.callback_PickPlotMarker,'Next Marker Type');
                in_addmenu(menus, 0 ,@gui.callback_PickColorMap,'Next Colormap');
                in_addmenu(menus, 1 ,@in_cleanumap,'Clean tSNE/UMAP/PHATE Plot');
            case 5
                menus = uimenu(FigureHandle, 'Text', 'Anno&tate');
                in_addmenu(menus, 0, {@in_DetermineCellTypeClustersGeneral, true}, "Annotate Cell Types Using PanglaoDB Marker Genes");
                in_addmenu(menus, 0, {@in_DetermineCellTypeClustersGeneral, false}, 'Annotate Cell Types Using Customized Marker Genes...');
                in_addmenu(menus, 1, {@in_MergeCellSubtypes, 1}, 'Import Cell Annotation from SCE in Workspace...');
                in_addmenu(menus, 0, {@in_MergeCellSubtypes, 2}, 'Import Cell Annotation from SCE Data File...');
                in_addmenu(menus, 1, @in_Brush4Celltypes, "Annotate Cell Types for Brushed Cells");
                in_addmenu(menus, 0, @gui.callback_Brush4Markers, "Find Marker Genes for Brushed Cells");
                in_addmenu(menus, 0, @gui.callback_FindAllMarkers, "Make Marker Gene Heatmap");
                in_addmenu(menus, 1, {@in_CellCyclePotency, 1}, 'Estimate Cell Cycle Phase...');
                in_addmenu(menus, 0, {@in_CellCyclePotency, 2}, 'Estimate Differentiation Potency...');
                in_addmenu(menus, 0, {@in_CellCyclePotency, 3}, 'Estimate Stemness...');
                in_addmenu(menus, 0, {@in_CellCyclePotency, 4}, 'Estimate Dissociation Gene Ratio...');
                in_addmenu(menus, 1, @in_SingleClickSolution, 'Single Click Solution (from Raw Data to Annotation)...');
            case 6
                menus = uimenu(FigureHandle, 'Text', '&Analyze');
                in_addmenu(menus, 0, @gui.callback_CalculateGeneStats, 'Gene Expression (Statistics) Analysis...');
                in_addmenu(menus, 0, @in_EnrichrHVGs, 'Gene Variability (HVG Function) Analysis...');
                in_addmenu(menus, 0, @in_CompareCellScoreBtwCls, 'Gene Program (Cell Score) Analysis...');
                % in_addmenu(menus, 0, @gui.callback_SGSAnalysis, 'Stochastic Gene Silencing (SGS) Analysis [PMID:39828795]...');
                in_addmenu(menus, 1, @gui.callback_DEGene2Groups, 'Differential Expression (DE) Analysis...','D');
                in_addmenu(menus, 0, @gui.callback_DVGene2Groups, 'Differential Variability (DV) Analysis...','V');
                in_addmenu(menus, 0, @gui.callback_DPGene2Groups, 'Differential Program (DP) Analysis...','P');
                in_addmenu(menus, 1, @gui.callback_DEGene2GroupsBatch, 'DE Analysis in Cell Type Batch Mode...');
                in_addmenu(menus, 0, @gui.callback_DVGene2GroupsBatch, 'DV Analysis in Cell Type Batch Mode...');
                in_addmenu(menus, 0, @gui.callback_DPGene2GroupsBatch, 'DP Analysis in Cell Type Batch Mode...');
                in_addmenu(menus, 1, @gui.callback_DEVP2GroupsBatch, 'All (DE, DV, DP) Analysis in Cell Type Batch Mode...');                
                % in_addmenu(menus, 1, @gui.callback_RunEnrichr, 'Enrichr Analysis...');                
            case 7
                menus = uimenu(FigureHandle, 'Text', '&Network');
                in_addmenu(menus, 0, @in_Select5000Genes, 'Remove Less Informative Genes to Reduce Gene Space...');
                in_addmenu(menus, 0, @gui.callback_DrawNetwork, 'Plot GRN from Edge (Gene Pair) List...');
                in_addmenu(menus, 1, @gui.callback_BuildGeneNetwork, 'Build GRN with Selected Genes...');
                in_addmenu(menus, 0, @gui.callback_CompareGeneNetwork, 'Build & Compare GRNs...');
                in_addmenu(menus, 1, {@in_scTenifoldNet,1}, 'Construct GRN with All Genes - scTenifoldNet [PMID:33336197] 🐢...');
                in_addmenu(menus, 0, {@in_scTenifoldNet,2}, 'Construct & Compare GRNs - scTenifoldNet [PMID:33336197] 🐢...');
                in_addmenu(menus, 1, @gui.callback_scTenifoldKnk1, 'Virtual Gene Knockout - scTenifoldKnk [PMID:35510185] 🐢 ...');
                in_addmenu(menus, 0, @gui.callback_VirtualKOGenKI, 'Virtual Gene Knockout (GenKI/🐍) [PMID:37246643] ...');
                %in_addmenu(menus, 0, @gui.callback_VirtualKOGenKI, 'Virtual Gene Knockout (GenKI/🐍) [PMID:37246643] 🐢 ...');
                %in_addmenu(menus, 1, @gui.callback_scTenifoldXct, 'Cell-Cell Communication (scTenifoldXct/🐍) [PMID:36787742] 🐢 ...');
                in_addmenu(menus, 1, @gui.callback_scTenifoldXct, 'One-Sample Cell-Cell Communication (scTenifoldXct/🐍) [PMID:36787742]...');
                in_addmenu(menus, 0, @gui.callback_scTenifoldXct2, 'Two-Sample Cell-Cell Communication (scTenifoldXct2/🐍) [PMID:36787742]...');
                in_addmenu(menus, 0, @gui.callback_scTenifoldCko, 'Virtual Cell-Cell Communication Knockout - scTenifoldCko/🐍 [Experimental] 🐢 ...');
            case 8
                menus = uimenu(FigureHandle, 'Text', 'E&xternal');
                in_addmenu(menus, 0, @gui.i_resetrngseed, 'Set Random Seed...');
                in_addmenu(menus, 0, @gui.i_setextwd, 'Set Working Folder...');
                in_addmenu(menus, 1, @gui.i_setscimilaritymodelpath, '🔢 - Set Scimilarity Model...');
                in_addmenu(menus, 0, @gui.i_setllmmodel, '🤖 - Set Large Language Model...');
                in_addmenu(menus, 1, @gui.callback_RunEnrichr, '🌐 - Enrichr Analysis...')
                in_addmenu(menus, 1, @gui.i_setrenv, 'Set Up R (ℝ) Environment');
                in_addmenu(menus, 0, @gui.i_setpyenv, 'Set Up Python (🐍) Environment');
                in_addmenu(menus, 1, @in_RunSeuratWorkflow, 'ℝ - Run Seurat Workflow (Seurat) [PMID:25867923]...');
                in_addmenu(menus, 0, @in_RunMonocle3, 'ℝ - Pseudotime Analysis (Monocle3) [PMID:28825705]...');
                in_addmenu(menus, 0, {@in_CellCyclePotency, 5}, 'ℝ - Aneuploid/Diploid Analysis (copykat) [PMID:33462507]...');
                in_addmenu(menus, 0, @in_DecontX, 'ℝ - Detect Ambient RNA Contamination (DecontX) [PMID:32138770]...');
                % in_addmenu(menus, 1, @in_RunDataMapPlot, '🐍 - Run DataMapPlot (datamapplot)...');
                in_addmenu(menus, 0, @gui.callback_RunMemento, '🐍 - Memento DE/DV Analysis [PMID:39454576]...');
                in_addmenu(menus, 0, @in_DoubletDetection, '🐍 - Detect Doublets (Scrublet) [PMID:30954476]...');
                in_addmenu(menus, 0, @in_HarmonyPy, '🐍 - Batch Integration (Harmony) [PMID:31740819]...');
                in_addmenu(menus, 0, @in_SCimilarity, '🐍 - Annotate Cell Types (Scimilarity) [PMID:39566551]...');
                in_addmenu(menus, 0, {@in_SubsampleCells, 2}, '🐍 - Geometric Sketching (geosketch) [PMID:31176620]...');
                in_addmenu(menus, 0, @gui.callback_MELDPerturbationScore, '🐍 - MELD Perturbation Score (MELD) [PMID:33558698]...');
                
                % in_addmenu(menus, 1, @gui.callback_ExploreCellularCrosstalk, 'Talklr Intercellular Crosstalk [DOI:10.1101/2020.02.01.930602]...');
            case 9
                %menus = uimenu(FigureHandle, 'Text', '&Options');
                %in_addmenu(menus, 0, @gui.i_resetrngseed, 'Set Random Seed...');
                %in_addmenu(menus, 0, @gui.i_setextwd, 'Set Working Folder...');
                %in_addmenu(menus, 1, @gui.i_setllmmodel, 'Set LLM Provider && Model...');
            case 10
                menus = uimenu(FigureHandle, 'Text', '&Help');
                in_addmenu(menus, 0, {@(~, ~) web('https://scgeatoolbox.readthedocs.io/en/latest/')}, 'Online Documentation...');
                % in_addmenu(menus, 0, {@(~, ~) gui.gui_showrefinfo('Quick Installation',FigureHandle)}, 'Quick Installation Guide...');
                in_addmenu(menus, 0, {@(~, ~) gui.gui_showrefinfo('Shortcuts Guide',FigureHandle)}, 'Shortcuts User Guide...');
                in_addmenu(menus, 1, {@(~, ~) web('https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox-single-cell-gene-expression-analysis-toolbox')}, 'View scGEAToolbox on File Exchange...');
                in_addmenu(menus, 0, {@(~, ~) web('https://pubmed.ncbi.nlm.nih.gov/31697351/')}, 'Cite scGEAToolbox Paper...');
                in_addmenu(menus, 0, {@(~, ~) web('https://scholar.google.com/scholar?cites=4661048952867744439&as_sdt=5,44&sciodt=0,44&hl=en')}, 'Papers Citing scGEAToolbox...');
                in_addmenu(menus, 1, {@(~, ~) web('https://scgeatool.github.io/')}, 'Visit SCGEATOOL-Standalone Website...');
                in_addmenu(menus, 0, {@(~, ~) web('https://matlab.mathworks.com/open/github/v1?repo=jamesjcai/scGEAToolbox&file=online_landing.m')}, 'Run SCGEATOOL in MATLAB Online...');
                in_addmenu(menus, 0, @gui.callback_InstallApp, 'Install scgeatoolApp...');
                in_addmenu(menus, 1, @gui.callback_CheckUpdates, 'Check for Updates...');
                [~, ~, ~, im] = i_majvercheck([false false false true]);
                in_addmenu(menus, 1, {@(~, ~) gui.sc_simpleabout(FigureHandle, im)}, 'About SCGEATOOL');
        end
    end

    function createPushButtons(FigureHandle)
        in_addbuttonpush(0, 0, [], [], "");
        in_addbuttonpush(0, 0, @gui.callback_MultiGroupingView, "visibility_18dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Multi-grouping View...");
        in_addbuttonpush(0, 0, @gui.callback_CrossTabulation, "full_stacked_bar_chart_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Cross tabulation");
        in_addbuttonpush(0, 0, @gui.callback_ShowGeneExprGroup, "dataset_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Gene expression in groups");
        in_addbuttonpush(0, 1, @gui.callback_Dotplot, "icon-mat-blur-linear-10.gif", "Gene Expression Dot Plot...");
        in_addbuttonpush(0, 0, @gui.callback_Heatmap, "icon-mat-apps-20.gif", "Gene Expression Heatmap...");
        in_addbuttonpush(0, 0, @gui.callback_ScatterStemPlot, "signal_cellular_pause_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Gene Expression/Cell State Stem Plot...");
        in_addbuttonpush(0, 0, @gui.callback_Violinplot, "violinplot.gif", "Gene Expression/Cell State Violin Plot...");
        in_addbuttonpush(0, 0, @gui.callback_ScatterCorrPlot, "icon-mat-blur-off-10a.gif", "Correlation Plot...");
        in_addbuttonpush(0, 0, [], [], "");
        in_addbuttonpush(0, 1, @in_CompareCellScoreBtwCls, "barcode_scanner_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Cell score analysis--obtaining gene signature score for each cell");
        in_addbuttonpush(0, 0, @gui.callback_GetCellSignatureMatrix, "hexagon_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Cell state analysis--obtaining multiple gene signature scores to reveal functional state of cells");
        in_addbuttonpush(0, 1, @gui.callback_DEGene2Groups, "graphic_eq_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Differential expression (DE) analysis");
        in_addbuttonpush(0, 0, @gui.callback_DVGene2Groups, "edit_audio_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Differential variability (DV) analysis");
        in_addbuttonpush(0, 0, @gui.callback_DPGene2Groups, "data_exploration_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Differential program (DP) analysis");
        in_addbuttonpush(0, 0, [], [], "");
        in_addbuttonpush(0, 1, @gui.callback_BuildGeneNetwork, "graph_3_18dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Build gene regulatory network");
        in_addbuttonpush(0, 0, @gui.callback_CompareGeneNetwork, "graph_5_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Compare two scGRNs");
        in_addbuttonpush(0, 1, {@gui.i_savemainfig, 3}, "presentation.jpg", 'Save Figure to PowerPoint File...');
        gui.gui_3dcamera(DeftToolbarHandle, '', false, FigureHandle);
        pt = uitoggletool(DeftToolbarHandle);
        try
            load(fullfile(mfolder, 'resources', 'Images', 'colorbarcdata.mat'), 'CData');
            pt.CData = CData;
        catch
            pt.CData = rand(16, 16, 3);
        end
        pt.Tooltip = 'Insert Colorbar';
        pt.ClickedCallback = @in_addcolorbar;
        pt.Tag = "figToglColorbar";        
        in_addbuttonpush(0, 0, @gui.callback_CloseAllOthers, "xmark-square.jpg", "Close all other figures");
        in_addbuttonpush(0, 0, {@gui.i_resizewin, FigureHandle}, 'scale-frame-reduce.jpg', 'Resize Plot Window')
        in_addbuttontoggle(1, 0, {@in_togglebtfun, @in_turnoffuserguiding, "icon-mat-unfold-more-10.gif", "icon-mat-unfold-less-10.gif", false, "Turn on/off user onboarding toolbar"});
        in_addbuttonpush(1, 0, @gui.callback_ShowGeneExpr, "google-docs.jpg", "Select genes to show expression")
        in_addbuttonpush(1, 0, @in_ShowCellStates, "bookmark-book.jpg", "Show cell states")
        in_addbuttonpush(1, 0, @in_SelectCellsByQC, "filter-alt.jpg", "Filter genes & cells")
        in_addbuttonpush(1, 1, @in_labelcellgroups, "label.jpg", "Label cell groups");
        in_addbuttonpush(1, 0, @in_Brushed2NewCluster, "substract.jpg", "Add brushed cells to a new group")
        in_addbuttonpush(1, 0, @in_Brushed2MergeClusters, "union.jpg", "Merge brushed cells to same group")
        in_addbuttonpush(1, 0, @in_RenameCellTypeBatchID, "edit.jpg", "Rename cell type or batch ID");
        in_addbuttonpush(1, 0, @in_SingleClickSolution, "apple-swift.jpg", "Single-click cell type annotation")
        in_addbuttonpush(1, 0, [], [], "");
        in_addbuttonpush(1, 1, @in_ClusterCellsS, "bubble_chart_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Clustering using cell embedding (S)")
        in_addbuttonpush(1, 0, {@in_DetermineCellTypeClustersGeneral, true}, "auto_transmission_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Assign cell types to groups")
        in_addbuttonpush(1, 0, @in_Brush4Celltypes, "color-picker.jpg", "Assign cell type to selected cells");
        in_addbuttonpush(1, 1, @gui.callback_Brush4Markers, "medal1st.jpg", "Marker genes of brushed cells");
        in_addbuttonpush(1, 0, @gui.callback_FindAllMarkers, "grain_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Marker gene heatmap");
        in_addbuttonpush(1, 0, [], [], "");
        in_addbuttonpush(1, 1, @gui.callback_ShowClustersPop, "view-grid.jpg", "Show cell clusters/groups individually");
        in_addbuttonpush(1, 0, @gui.callback_SelectCellsByClass, "checklist_rtl_18dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Select cells by class/group & open in new window");
        in_addbuttonpush(1, 0, @in_DeleteSelectedCells, "variable_remove_18dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Delete brushed/selected cells");
        in_addbuttonpush(1, 0, @gui.callback_SaveX, "floppy-disk.jpg", "Export & save data");
        in_addbuttonpush(1, 1, @in_EmbeddingAgain, "cinema-old.jpg", "Embedding (tSNE, UMP, PHATE)");
        in_addbuttonpush(1, 0, @in_Switch2D3D, "perspective-view.jpg", "Switch 2D/3D");
        in_addbuttonpush(1, 0, @gui.callback_PickPlotMarker, "palette.jpg", "Switch scatter plot marker type");
        in_addbuttonpush(1, 0, @gui.callback_PickColorMap, "color-wheel.jpg", "Pick new color map");
        in_addbuttonpush(1, 0, @in_RefreshAll, "refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", "Refresh");
        in_addbuttonpush(2, 0, @in_turnonuserguiding, "icon-fa-thumb-tack-10.gif", "Turn on user guiding toolbar");
        in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_SelectCellsByQC, "icon-mat-filter-1-10.gif", "filter-alt.jpg", true, "Filter genes & cells"});
        in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_EmbeddingAgain, "icon-mat-filter-2-10.gif", "cinema-old.jpg", true, "Embedding (tSNE, UMP, PHATE)"});
        in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_ClusterCellsS, "icon-mat-filter-3-10.gif", "bubble_chart_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", true, "Clustering using embedding S"});
        in_addbuttontoggle(2, 0, {@in_togglebtfun, @in_DetermineCellTypeClustersGeneral, "icon-mat-filter-4-10.gif", "auto_transmission_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", true, "Assign cell types to groups"});
        in_addbuttontoggle(2, 0, {@in_togglebtfun, @gui.callback_SaveX, "icon-mat-filter-5-10.gif", "floppy-disk.jpg", true, "Export & save data"});
    end

    function in_sc_openscedlg(~, event)        
        if strcmp(event.EventName,'KeyPress') && ~ismember(event.Key,{'return','space','i','I'}), return; end
        clickType = get(FigureHandle, 'SelectionType');
        if strcmp(clickType,'alt'), return; end
        
        set(button1,'Enable','off');
        if ~isempty(sce) && sce.NumCells > 0            
            if ~strcmp(gui.myQuestdlg(FigureHandle,'Current SCE will be replaced. Continue?',''),'Yes') 
                set(button1,'Enable','on');
                return;
            end
        end
        [sce] = gui.sc_openscedlg([],[],FigureHandle);
        if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end
        if ~isempty(sce) && sce.NumCells > 0 && sce.NumGenes > 0
            guidata(FigureHandle, sce);
            c=[];
            in_RefreshAll([], [], false, false);
        else
            set(button1,'Enable','on');
            if ~gui.i_isuifig(FigureHandle)
                uicontrol(button1);
            end
            if ~isempty(sce)
                gui.myWarndlg(FigureHandle, 'Imported SCE contains no cells.');
            end
        end
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

    function in_cleanumap(~, ~)
        answer = gui.myQuestdlg(FigureHandle,'Select embedding method label.', ...
            '',{'tSNE','UMAP','PHATE'},'tSNE');
        if isempty(answer), return; end
        a = colormap;
        gui.i_baredrplot(hAx, [], answer, FigureHandle);
        colormap(hAx,a);
    end

    function in_CompareCellScoreBtwCls(src, events)
        if gui.callback_CompareCellScoreBtwCls(src, events)
           sce = guidata(FigureHandle);
        end
    end

    function in_CellCyclePotency(src, events, typeid)
        if gui.callback_CellCyclePotency(src, events, typeid)
           sce = guidata(FigureHandle);
        end
    end

    function in_RunMonocle3(src, events)
        if gui.callback_RunMonocle3(src, events)
            sce = guidata(FigureHandle);
            [y, idx] = ismember({'monocle3_pseudotime'}, ...
                 sce.list_cell_attributes(1:2:end));
            if y
                answer = gui.myQuestdlg(FigureHandle,'Color cells using pseudotime T and show Monocle3 embedding S?',...
                    '',{'Yes','Color cells only','Show embedding only'},'Yes');
                switch answer
                    case 'Yes'
                        sce.c = sce.list_cell_attributes{idx*2};
                        sce.s = sce.struct_cell_embeddings.('monocle2d');
                        [c, cL] = grp2idx(sce.c);
                        in_RefreshAll(src, [], true, false);
                    case 'Color cells only'
                        sce.c = sce.list_cell_attributes(idx*2);
                        [c, cL] = grp2idx(sce.c);
                        in_RefreshAll(src, [], true, false);
                    case 'Show embedding only'
                        sce.s = sce.struct_cell_embeddings.('monocle2d');
                        in_RefreshAll(src, [], true, false);
                    otherwise
                end                
            end
        end
    end


    function in_turnonuserguiding(~, ~)
        % setpref('scgeatoolbox','useronboardingtoolbar',true);
        % set(UserToolbarHandle, 'Visible', 'on');
        Button = gui.gui_userguidingpref(false);
        switch Button
            case 'Yes'
            case 'No'
                in_turnoffuserguiding;
            case 'Cancel'
        end
    end

    function in_turnoffuserguiding(~, ~)
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
            %answer=gui.myQuestdlg(FigureHandle, 'Show User Onboarding Toolbar again next time?','');
            %switch answer
            %    case 'Yes'
            %        setpref('scgeatoolbox','useronboardingtoolbar',true);
            %    case 'No'
            %        setpref('scgeatoolbox','useronboardingtoolbar',false);
            %end
        end
    end

    function in_addmenu(menuHdl, sepTag, callbackFnc, tooltipTxt, acchar)
        if nargin<5, acchar=''; end
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
            'Callback', callbackFnc, ...
            'Accelerator', acchar, ...
            'Tag', "figMenu" + matlab.lang.makeValidName(tooltipTxt));
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
       
        if numel(ptImgCell) >= ptimgidx && ~isempty(ptImgCell{1, ptimgidx})
            pkg.i_addbutton2fig(barhandle, sepTag, ...
                callbackFnc, ptImgCell{1, ptimgidx}, tooltipTxt);
        else
            [~, ptImgCell{1, ptimgidx}] = pkg.i_addbutton2fig(barhandle, ...
                sepTag, callbackFnc, imgFil, tooltipTxt);
            ptImgCell{2, ptimgidx} = tooltipTxt;
        end
        ptimgidx = ptimgidx + 1;
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
        %callbackFnc =
        %    1×6 cell array
        %    {function_handle}    {function_handle}    {["icon-mat-filt…"]}    {["plotpicker-co…"]}    {[1]}    {["Assign cell t…"]}
        %pt.Tag = "figTogl" + matlab.lang.makeValidName(tooltipTxt);
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
                try
                    gui.myHelpdlg(FigureHandle, sprintf('%s\n%s', upper(tooltipTxt), s));
                catch
                end
            end
        else
            func(src);
        end
    end

    function [ptImage] = in_getPtImage(imgFil)
        try
            [ptImage, map] = imread(fullfile(mfolder, 'resources', 'Images', imgFil));
        catch
            try
                [ptImage, map] = imread(fullfile(matlabroot,'toolbox', ...
                    'matlab','icons', imgFil));
            catch
                map = [];
                ptImage = rand(16, 16, 3); 
            end
        end
        if ~isempty(map), ptImage = ind2rgb(ptImage, map); end
        if size(ptImage, 3) == 1, ptImage = cat(3, ptImage, ptImage, ptImage); end
    end

% ------------------------
% Callback Functions
% ------------------------

    function in_addcolorbar(~,~)
        cbtogg = findall(FigureHandle, 'Tag', 'figToglColorbar');
        if ~isempty(cbtogg) && isequal(cbtogg,gcbo) && strcmpi(get(cbtogg,'State'),'on')
            colorbar(hAx);
        else
            colorbar(hAx,'off')
        end
    end

    function in_closeRequest(hObject, ~)
        if ~(ismcc || isdeployed)
            if isempty(sce) || sce.NumCells==0
                ButtonName = 'no';
            else
                ButtonName = gui.myQuestdlg(FigureHandle, 'Save SCE before closing SCGEATOOL?');
            end
            switch lower(ButtonName)
                case 'yes'
                    if ~isempty(callinghandle)
                        guidata(callinghandle, sce);
                        delete(hObject);
                        gui.myHelpdlg(FigureHandle, 'SCE updated.');
                    else
                        if gui.callback_SaveX(FigureHandle,[])
                            pause(1);
                            delete(hObject);
                        end
                    end
                case 'cancel'
                    return;
                case 'no'
                    delete(FigureHandle);
                otherwise
                    return;
            end
        else
            delete(hObject);
        end
    end

    function in_qcviolin(~, ~)
        gui.i_qcviolin(sce.X, sce.g, FigureHandle);
    end

    function in_RunDataMapPlot(src, ~)
        if ~pkg.i_checkpython
            gui.myWarndlg(FigureHandle, 'Python is not installed.');
            return;
        end
        ndim = 2;
        [vslist] = gui.i_checkexistingembed(sce, ndim);
        if isempty(h.ZData) && size(sce.s,2)==2 && length(vslist) <= 1
            gui.callback_RunDataMapPlot(src, []);
        elseif isempty(h.ZData) && size(sce.s,2)==2 && length(vslist) > 1            
            switch gui.myQuestdlg(FigureHandle, 'Using current 2D embedding?')
                case 'Yes'
                    gui.callback_RunDataMapPlot(src, []);
                case 'No'
                    [sx] = gui.i_pickembedvalues(sce, 2, FigureHandle);
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
            if strcmp(gui.myQuestdlg(FigureHandle, 'This function requires 2D embedding. Continue?'), 'Yes'), in_Switch2D3D(src,[]); end
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
            switch gui.myQuestdlg(FigureHandle, 'Import annotation for all cells or just cells of a subtype?', '', ...
                {'All Cells', 'Subtype Cells', 'Cancel'}, 'All Cells')
                case 'All Cells'
                    allcell = true;
                case 'Subtype Cells'
                    allcell = false;
                case 'Cancel'
                    return;
                otherwise
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
            if sce.NumCells==0
                gui.myWarndlg(FigureHandle, 'Merged SCE contains no cells.');
                return;
            else
                in_RefreshAll(src, [], true, false);
                gui.myHelpdlg(FigureHandle, sprintf('%s SCEs merged.', upper(s)));
            end
        end
    end

    function in_WorkonSelectedGenes(src, ~, type)
        if nargin < 1, type = 'name'; end

        switch type
            case 'name'
                [glist] = gui.i_selectngenes(sce, [], FigureHandle);
                if isempty(glist), return; end
                [y, idx] = ismember(glist, sce.g);
                if ~all(y)
                    gui.myErrordlg(FigureHandle, 'Runtime error.','');
                    return;
                end
            case 'hvg'
                k = gui.i_inputnumk(2000, 1, sce.NumGenes, 'the number of HVGs');
                if isempty(k), return; end
                answer = gui.myQuestdlg(FigureHandle, 'Which HVG detecting method to use?', '', ...
                    {'Splinefit Method [PMID:31697351]', ...
                    'Brennecke et al. (2013) [PMID:24056876]'}, ...
                    'Splinefit Method [PMID:31697351]');                
                switch answer
                    case 'Splinefit Method [PMID:31697351]'
                        fw = gui.myWaitbar(FigureHandle);
                        T = sc_splinefit(sce.X, sce.g);
                    case 'Brennecke et al. (2013) [PMID:24056876]'
                        fw = gui.myWaitbar(FigureHandle);
                        T = sc_hvg(sce.X, sce.g);                        
                    otherwise
                        return;
                end
                glist = T.genes(1:min([k, sce.NumGenes]));
                [y, idx] = ismember(glist, sce.g);
                if ~all(y)
                    gui.myErrordlg(FigureHandle, 'Runtime error.','');
                    return;
                end
            case 'ligandreceptor'
                fw = gui.myWaitbar(FigureHandle);
                load(fullfile(mfolder, 'resources', 'Ligand_Receptor', ...
                     'Ligand_Receptor_more.mat'), 'ligand','receptor');
                idx = ismember(upper(sce.g), unique([ligand; receptor]));
                if ~any(idx)
                    gui.myErrordlg(FigureHandle, 'Runtime error: No gene left after selection.','');
                    return;
                end
                if sum(idx) < 50                    
                    if ~strcmp(gui.myQuestdlg(FigureHandle, 'Few genes (n < 50) selected. Continue?',''), 'Yes'), return; end
                end
        end

        sce.g = sce.g(idx);
        sce.X = sce.X(idx, :);
        gui.myWaitbar(FigureHandle, fw);
        guidata(FigureHandle, sce);
        in_RefreshAll(src, [], true, false);
    end

    function in_SubsampleCells(src, ~, methodoption)
        if nargin < 3, methodoption = []; end        
        if ~strcmp(gui.myQuestdlg(FigureHandle, 'This function subsamples 50% of cells. Continue?',''), 'Yes'), return; end
        if isempty(methodoption)
            answer = gui.myQuestdlg(FigureHandle, 'Select method:', '', ...
                {'Uniform Sampling', ...
                'Geometric Sketching [PMID:31176620]'}, 'Uniform Sampling');
            switch answer
                case 'Uniform Sampling'
                    methodoption = 1;
                case 'Geometric Sketching [PMID:31176620]'
                    if ~pkg.i_checkpython
                        gui.myWarndlg(FigureHandle, 'Python not installed.','');
                        return;
                    end
                    methodoption = 2;
                otherwise
                    return;
            end
        end

        tn = round(sce.NumCells/2);
        if methodoption == 1
            rng("shuffle");
            idx = randperm(sce.NumCells);
            ids = idx(1:tn);
        elseif methodoption == 2
            if ~gui.gui_showrefinfo('Geometric Sketching [PMID:31176620]', ...
                    FigureHandle), return; end
            fw = gui.myWaitbar(FigureHandle);
            Xn = log1p(sc_norm(sce.X))';
            [~, Xn] = pca(Xn, 'NumComponents', 300);
            gui.myWaitbar(FigureHandle, fw);
            try
                ids = run.py_geosketch(Xn, tn);
            catch ME
                gui.myWaitbar(FigureHandle, fw, true);
                gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
                return;
            end

        end
        if ~isempty(ids)
            sce = sce.selectcells(ids);
            c = sce.c;
            in_RefreshAll(src, [], true, false);
            guidata(FigureHandle, sce);
        else
            gui.myErrordlg(FigureHandle, 'Runtime error. No action is taken.');
        end
    end

    function in_SingleClickSolution(src, ~)
        if ~isprop(sce, 'c_cell_type_tx')
            disp('The sce object does not have the property c_cell_type_tx.');
            return;
        end
        if ~all(sce.c_cell_type_tx == "undetermined")            
            if ~strcmp(gui.myQuestdlg(FigureHandle, ...
                    "Your data has been embedded and annotated. " + ...
                    "Single Click Solution will re-embed and " + ...
                    "annotate cells. Current embedding and " + ...
                    "annotation will be overwritten. Continue?"), 'Yes')
                return; 
            end
        else
            if ~gui.gui_showrefinfo('Single Click Solution', ...
                    FigureHandle), return; end
        end
        %if isempty(speciestag)
            speciestag = gui.i_selectspecies(2, false, FigureHandle);
        %end
        if isempty(speciestag), return; end

        prompt = {
            'tSNE Embedding?', ...            
            'Add UMAP Embedding?', ...            
            'Add PHATE Embedding?', ...
            'Estimate Cell Cycles?', ...
            'Estimate Differentiation Potency of Cells?'};
        dlgtitle = '';
        dims = [1, 85];        
        definput = {'Yes', 'No', 'No', 'Yes', 'No'};

        if gui.i_isuifig(FigureHandle)
            answer = gui.myInputdlg(prompt, dlgtitle, definput, FigureHandle);
        else
            answer = inputdlg(prompt, dlgtitle, dims, definput);
        end

        if isempty(answer)
            return;
        end

        fw = gui.myWaitbar(FigureHandle);
        gui.myWaitbar(FigureHandle, fw, false, '','Basic QC Filtering...', 1/8);
        sce = sce.qcfilter;

        count = 1;
        if ~(strcmpi(answer{count},'Yes') || strcmpi(answer{count},'Y'))
            gui.myErrordlg(FigureHandle, 'tSNE Embedding has to be included.','');
            return;        
        end

        count = count + 1;
        if strcmpi(answer{count},'Yes') || strcmpi(answer{count},'Y')
            gui.myWaitbar(FigureHandle, fw, false, '', 'Embeding Cells Using UMAP...', 2/8);
            sce = sce.embedcells('umap3d', true, true, 3);
        end
        count = count + 1;
        if strcmpi(answer{count},'Yes') || strcmpi(answer{count},'Y')
            gui.myWaitbar(FigureHandle, fw, false, '', 'Embeding Cells Using PHATE...', 2/8);
            sce = sce.embedcells('phate3d', true, true, 3);
        end
        gui.myWaitbar(FigureHandle, fw, false, '', 'Embeding Cells Using tSNE...', 2/8);
        try
            sce = sce.embedcells('tsne3d', true, true, 3);
        catch ME
            gui.myWaitbar(FigureHandle, fw);
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        gui.myWaitbar(FigureHandle, fw, false, '', 'Clustering Cells Using K-means...', 3/8);
        
        sce = sce.clustercells([], [], true);
        gui.myWaitbar(FigureHandle, fw, false, '', 'Annotating Cell Types Using PanglaoDB...', 4/8);
        sce = sce.assigncelltype(speciestag, false);

        count = count + 1;
        if strcmpi(answer{count},'Yes') || strcmpi(answer{count},'Y')
            gui.myWaitbar(FigureHandle, fw, false, '', 'Estimate Cell Cycles...', 5/8);
            sce = sce.estimatecellcycle;
        end

        count = count + 1;
        if strcmpi(answer{count},'Yes') || strcmpi(answer{count},'Y')
            gui.myWaitbar(FigureHandle, fw, false, '', ...
                'Estimate Differentiation Potency of Cells...', 6/8);
            sce = sce.estimatepotency(speciestag);
        end
        gui.myWaitbar(FigureHandle, fw, false, '', '', 7/8);
        [c,cL] = grp2idx(sce.c_cell_type_tx);
        sce.c = c;
        gui.myWaitbar(FigureHandle, fw);
        in_RefreshAll(src, [], true, false);
        ix_labelclusters(true);
        %setappdata(FigureHandle, 'cL', cL);
        guidata(FigureHandle, sce);
    end

    function in_SelectCellsByQC(src, ~)
        %oldsce = sce;
        % oldn = sce.NumCells;
        % oldm = sce.NumGenes;
        sce.c = c;
        guidata(FigureHandle, sce);
        try
            [requirerefresh, highlightindex] = ...
                gui.callback_SelectCellsByQC(src);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        if requirerefresh
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c);
            in_RefreshAll(src, [], true, false);
            % newn = sce.NumCells;
            % newm = sce.NumGenes;
            % answer = gui.myQuestdlg(FigureHandle, sprintf('%d cells removed; %d genes removed.', ...
            %         oldn-newn, oldm-newm),'','Accept Changes', 'Undo Changes', 'Accept Changes');
            % if ~strcmp(answer, 'Accept Changes')
            %     sce = oldsce;
            %     [c, cL] = grp2idx(sce.c);
            %     in_RefreshAll(src, [], true, false);
            %     guidata(FigureHandle, sce);
            % end
        end
        if ~isempty(highlightindex)
            h.BrushData = highlightindex;
        end
    end

    function in_Select5000Genes(src, ~)
        oldm = sce.NumGenes;
        oldn = sce.NumCells;
        [requirerefresh, scenew] = gui.callback_Select5000Genes(src);
        if requirerefresh
            sce = scenew;
            c=sce.c;
            
            newm = sce.NumGenes;
            newn = sce.NumCells;
            gui.myHelpdlg(FigureHandle, sprintf('%d cells removed; %d genes removed.', ...
                oldn-newn, oldm-newm));
            guidata(FigureHandle, sce);   
            in_RefreshAll(src, [], true, false);
        end
    end

    function in_RunSeuratWorkflow(src, ~)
        extprogname = 'R_Seurat';
        preftagname = 'externalwrkpath';
        [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
        if isempty(wkdir), return; end
        [ok] = gui.i_confirmscript('Run Seurat Workflow (Seurat/R)?', ...
            'R_Seurat', 'r');
        if ~ok, return; end

        [ndim] = gui.i_choose2d3d;
        if isempty(ndim), return; end
        fw = gui.myWaitbar(FigureHandle);
        try
            [sce] = run.r_seurat(sce, ndim, wkdir, true);
            [c, cL] = grp2idx(sce.c);
        catch
            gui.myWaitbar(FigureHandle, fw);
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
        guidata(FigureHandle, sce);
        in_RefreshAll(src, [], true, false);
    end

    function in_DecontX(~, ~)
        if ~gui.gui_showrefinfo('DecontX [PMID:32138770]', FigureHandle), return; end
        extprogname = 'R_decontX';
        preftagname = 'externalwrkpath';
        [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
        if isempty(wkdir), return; end
        fw = gui.myWaitbar(FigureHandle);
        try
            [Xdecon, ~] = run.r_decontX(sce, wkdir);
        catch ME
           gui.myWaitbar(FigureHandle, fw);
           gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
           return;
        end
        gui.myWaitbar(FigureHandle, fw);

        if ~(size(Xdecon, 1) == sce.NumGenes && size(Xdecon, 2) == sce.NumCells)
           gui.myErrordlg(FigureHandle, 'DecontX runtime error.','');
           return;
        end

%       figure('WindowStyle', 'modal');
%       gui.i_stemscatter(sce.s, contamination);
%       zlabel('Contamination rate')
%       title('Ambient RNA contamination')        
        if strcmp(gui.myQuestdlg(FigureHandle, "Remove contamination? Click ''No'' to keep data unchanged."), 'Yes')
            sce.X = round(Xdecon);
            guidata(FigureHandle, sce);
            gui.myHelpdlg(FigureHandle, 'Contamination removed.');
        end
    end

    function in_SCimilarity(src, events)
        if ~pkg.i_checkpython
            gui.myWarndlg(FigureHandle, 'Python not installed.');
            return;
        end
        if ~gui.gui_showrefinfo('SCimilarity [PMID:39566551]', ...
                FigureHandle), return; end
        if gui.callback_RunSCimilarity(src, events)
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c_cell_type_tx);
            in_RefreshAll(src, [], true, false);            
        end
    end

    function in_HarmonyPy(src, ~)
        if ~pkg.i_checkpython
            gui.myWarndlg(FigureHandle, 'Python is not installed.');
            return;
        end        
        if ~gui.gui_showrefinfo('Harmony [PMID:31740819]', FigureHandle), return; end
        if numel(unique(sce.c_batch_id)) < 2
            gui.myWarndlg(FigureHandle, 'No batch effect (SCE.C_BATCH_ID is empty)');
            return;
        end
        [c1] = grp2idx(sce.c);
        [c2] = grp2idx(sce.c_batch_id);
        if ~isequal(c1, c2)
            answer = gui.myQuestdlg(FigureHandle, 'Color cells by batch id (SCE.C_BATCH_ID)?', '');
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

            ButtonName = gui.myQuestdlg(FigureHandle, 'Update Saved Embedding?', '');
            switch ButtonName
                case 'Yes'
                    [methodtag] = gui.i_pickembedmethod(FigureHandle);
                    if isempty(methodtag), return; end
                    [ndim] = gui.i_choose2d3d;
                    if isempty(ndim), return; end
                    methoddimtag = sprintf('%s%dd',methodtag, ndim);
                    if ismember(methoddimtag, fieldnames(sce.struct_cell_embeddings))
                        sce.struct_cell_embeddings.(methodtag) = sce.s;
                    end
                    gui.myHelpdlg(FigureHandle, ...
                        sprintf('%s Embedding is updated.', methoddimtag));
            end
        end
        guidata(FigureHandle, sce);
    end

    function in_DoubletDetection(src, ~)
        if ~pkg.i_checkpython
            gui.myWarndlg(FigureHandle, 'Python not installed.');
            return;
        end
        if ~gui.gui_showrefinfo('Scrublet [PMID:30954476]', FigureHandle), return; end
        if numel(unique(sce.c_batch_id)) > 1            
            if ~strcmp(gui.myQuestdlg(FigureHandle, '"When working with data from multiple samples, run Scrublet on each sample separately." Your data contains multiple samples (cells with different c_batch_id). Continue?',''), 'Yes'), return; end
        end
        [isDoublet, doubletscore, methodtag, done] = gui.callback_DoubletDetection(src);
        if done && ~any(isDoublet)
            gui.myHelpdlg(FigureHandle, 'No doublet detected.');
            return;
        end
        if done && any(isDoublet) && sce.NumCells == length(doubletscore)
            tmpf_doubletdetection = figure('WindowStyle', 'modal');
            gui.i_stemscatter(sce.s, doubletscore);
            zlabel('Doublet Score')
            title(sprintf('Doublet Detection (%s)', methodtag))
            if strcmp(gui.myQuestdlg(FigureHandle, sprintf("Remove %d doublets?", sum(isDoublet))),'Yes')
                close(tmpf_doubletdetection);
                % i_deletecells(isDoublet);
                sce = sce.removecells(isDoublet);
                guidata(FigureHandle, sce);
                [c, cL] = grp2idx(sce.c);
                in_RefreshAll(src, [], true, false);
                gui.myHelpdlg(FigureHandle, 'Doublets deleted.');
            end
        end
    end

    function in_MergeSubCellTypes(src, ~)
        if isempty(sce.c_cell_type_tx), return; end
        newtx = erase(sce.c_cell_type_tx, "_{"+digitsPattern+"}");
        if isequal(sce.c_cell_type_tx, newtx)
            gui.myHelpdlg(FigureHandle, "No sub-clusters are meraged.");
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
        %if strcmpi(entag,'on')
            %set(tb1,'Visible','off');
            %set(tb2,'Visible','off');
            
            %px=FigureHandle.Position;
            %px(4)=px(4)-50;
            %FigureHandle.Position=px;
        %end
        set(DeftToolbarHandle,'Visible',entag);
        set(MainToolbarHandle,'Visible',entag);
        showuseronboarding = getpref('scgeatoolbox', ...
            'useronboardingtoolbar',false);
        switch entag
            case 'on'
                if showuseronboarding
                    set(UserToolbarHandle,'Visible','on');
                else
                    set(UserToolbarHandle,'Visible','off');
                end
                if isvalid(button1), set(button1,'Visible','off'); end
                if isvalid(button2), set(button2,'Visible','off'); end
            case 'off'
                set(UserToolbarHandle,'Visible','off');
        end
        menusv={m_file, m_edit, m_view, m_plot, m_anno, m_tool, m_ntwk,...
            m_extn};
        for j=1:length(menusv)
            a=allchild(menusv{j});
            for k=1:length(a)
                a(k).Enable=entag;
            end
        end
        a=allchild(m_file);
        a(end).Enable='on';
        a(end-1).Enable='on';        
        a=allchild(m_extn);
        a(end).Enable='on';
        a(end-1).Enable='on';
        a(end-2).Enable='on';
        a(end-3).Enable='on';
        a(end-4).Enable='on';
        a(end-5).Enable='on';
        a(end-6).Enable='on';
        %a=allchild(m_optn);
        %a(end).Enable='on';
        %a(end-1).Enable='on';
        %a(end-2).Enable='on';
    end

    function in_RefreshAll(src, ~, keepview, keepcolr)
        if nargin < 4, keepcolr = false; end
        if nargin < 3, keepview = false; end

        if gui.i_isuifig(FigureHandle)
            keepcolr = false;
            keepview = false;
        end

        if keepview || keepcolr
            [para] = gui.i_getoldsettings(src);
        end
        if isempty(sce) || sce.NumCells==0
            if ~isempty(h) && isvalid(h)
                delete(h);
            end
            set(hAx,'Visible','off');
            in_EnDisableMenu('off');
            return;
        end
        in_EnDisableMenu('on');
        %        if isvalid(button1), set(button1,'Visible','off'); end
        %        if isvalid(button2), set(button2,'Visible','off'); end
        
        if ~gui.i_isuifig(FigureHandle)
            if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end
        end
        
        % [c,cL]=grp2idx(sce.c);
        % was3d = ~isempty(h.ZData);        

        if isempty(c), [c,cL] = grp2idx(sce.c); end
        if size(sce.s, 2) >= 3
            if keepview, [ax, bx] = view(hAx); end
            h = gui.i_gscatter3(sce.s, c, methodid, [], hAx);
            %h=scatter3(hAx,rand(100,1),rand(100,1),rand(100,1));
            if keepview && bx~=90
                view(hAx, ax, bx);
            end
        else        % otherwise going to show 2D
            if keepview, [ax, bx] = view(hAx); end
            h = gui.i_gscatter3(sce.s(:, 1:2), c, methodid, [], hAx);
            if keepview && bx==90
                view(hAx, ax, bx);
            end
        end
        if keepview
            if isfield(para,'oldMarker')
                h.Marker = para.oldMarker;
            end
            if isfield(para,'oldSizeData')
                h.SizeData = para.oldSizeData;
            end
        end
        if keepcolr
            if isfield(para, 'oldColorMap')
                colormap(hAx,para.oldColorMap);
            end
        else
            kc = numel(unique(sce.c));
            colormap(hAx,pkg.i_mycolorlines(kc));
        end
        title(hAx, sce.title);
        subtitle(hAx, '[genes x cells]');
        hAx.Toolbar.Visible = 'on';
        set(hAx,'Visible','on');
    end

    function in_Switch2D3D(src, ~)
        [para] = gui.i_getoldsettings(hAx);
        [ax, bx] = view(hAx);        
        if bx == 90   % isempty(h.ZData)               % current 2D            
            if ~strcmp(gui.myQuestdlg(FigureHandle, 'Switch to 3D?',''), 'Yes'), return; end
            if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end
            if size(sce.s, 2) >= 3
                h = gui.i_gscatter3(sce.s, c, methodid, [], hAx);
                if ~isempty(ax) && ~isempty(bx) && ~any([ax, bx] == 0)
                    view(hAx, ax, bx);
                else
                    view(hAx, 3);
                end
            else
                [vslist] = gui.i_checkexistingembed(sce, 3);
                if isempty(vslist)
                    in_EmbeddingAgain(src, [], 3);
                else
                    ansx = gui.myQuestdlg(FigureHandle, 'Using existing 3D embedding? Select "No" to re-embed.');
                    switch ansx
                        case 'Yes'
                            [sx] = gui.i_pickembedvalues(sce, 3, FigureHandle);
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
            if ~strcmp(gui.myQuestdlg(FigureHandle, 'Switch to 2D?',''), 'Yes'), return; end
            [vslist] = gui.i_checkexistingembed(sce, 2);
            
            if ~isempty(vslist)
                answer = gui.myQuestdlg(FigureHandle, 'How to make 2D embedding?','', ...
                    {'Pick existing 2D','Re-embed cells', ...
                    'Reduce current 3D'},'Pick existing 2D');
            else
                answer = gui.myQuestdlg(FigureHandle, 'How to make 2D embedding?','', ...
                    {'Embed Cells to 2D', ...
                    'Project Current 3D Embedding to 2D','Cancel'},'Embed Cells to 2D');
            end

            switch answer
                case 'Cancel'
                    return;
                case 'Re-embed cells'
                    in_EmbeddingAgain(src, [], 2);
                case {'Project Current 3D Embedding to 2D', 'Reduce current 3D'}
                    [ax, bx] = view(hAx);
                    answer2 = gui.myQuestdlg(FigureHandle, 'Which view to be used to project 3D to 2D?', '', ...
                        {'X-Y Plane', 'Screen/Camera', 'PCA-Rotated'}, 'X-Y Plane');
                    switch answer2
                        case 'X-Y Plane'
                            sx = sce.s;
                        case 'Screen/Camera'
                            sx = pkg.i_3d2d(sce.s, ax, bx);
                        case {'PCA-Rotated'}
                            [~, sx] = pca(sce.s);
                        otherwise
                            return;
                    end
                    h = gui.i_gscatter3(sx(:, 1:2), c, ...
                        methodid, [], hAx);
                    sce.s = sx(:, 1:2);
                    title(hAx, sce.title);
                    subtitle(hAx, '[genes x cells]');
                    h.Marker = para.oldMarker;
                    h.SizeData = para.oldSizeData;
                    colormap(hAx,para.oldColorMap);
                    view(hAx, 2);
                    return;
                case 'Pick existing 2D'
                    [sx] = gui.i_pickembedvalues(sce, 2, FigureHandle);
                    if ~isempty(sx) && size(sx,1) == sce.NumCells
                        sce.s = sx;
                    else                        
                        return;
                    end
            end
        end
        guidata(FigureHandle, sce);
        
        in_RefreshAll(src, [], true, true);   % keepview, keepcolr
    end

    function in_AddEditCellAttribs(~,~)        
        switch gui.myQuestdlg(FigureHandle, 'Add or edit cell attribute?','',{'Add','Edit','Cancel'},'Add')
            case 'Edit'
                addnew = false;
            case 'Add'
                addnew = true;
            otherwise
                return;
        end
        [sceback, needupdate] = gui.sc_cellattribeditor(sce, addnew, FigureHandle);
        if needupdate
            sce = sceback;
            guidata(FigureHandle, sce);
        end
    end

    function in_ExportCellAttribTable(~,~)        
        T = pkg.makeattributestable(sce);
        gui.i_exporttable(T,true, ...
            "Tcellattrib","CellAttribTable",[],[],FigureHandle);
    end

    function in_RenameCellTypeBatchID(src, ~, answer)
        if nargin < 3 || isempty(answer)
            answer = gui.myQuestdlg(FigureHandle, 'Rename cell type, batch ID, or gene name?', ...
                '',{'Cell type', 'Batch ID', 'Others...'}, 'Cell type');
        end
        switch answer
            case 'Cell type'
                [requirerefresh] = gui.callback_RenameCellType(src);
            case 'Batch ID'
                [requirerefresh] = gui.callback_RenameBatchID(src);
            case 'Gene name'
                [requirerefresh] = gui.callback_RenameGenes(src);
            case 'Others...'
                [requirerefresh, answer] = gui.callback_RenameOthers(src);
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
                case 'Cluster ID'
                    [c, cL] = grp2idx(sce.c_cluster_id);
                otherwise
                    return;
            end
            ix_labelclusters(false);
        end
    end

    function in_EmbeddingAgain(src, ~, ndim)
        readytoshow = false;
        if nargin<3, ndim=[]; end
        if isempty(sce.struct_cell_embeddings)
            sce.struct_cell_embeddings = pkg.e_makeembedstruct;
        end
        if ~isfield(sce.struct_cell_embeddings, 'metaviz3d')
            sce.struct_cell_embeddings.('metaviz3d') = [];
        end
        [vslist] = gui.i_checkexistingembed(sce, ndim);
        if ~isempty(vslist)
            answer = gui.myQuestdlg(FigureHandle, 'Using exsiting embedding?');
        else
            answer = 'No';
        end
        if strcmp(answer, 'Yes')
            [sx] = gui.i_pickembedvalues(sce,[],FigureHandle);
            if ~isempty(sx), sce.s = sx; end
        elseif strcmp(answer, 'No')
            [methodtag] = gui.i_pickembedmethod(FigureHandle);
            if isempty(methodtag), return; end
            %if isempty(ndim), [ndim] = gui.i_choose2d3dnmore; end
            %if isempty(ndim), return; end
            %methoddimtag = sprintf('%s%dd',methodtag, ndim); 

            switch gui.myQuestdlg(FigureHandle, 'Overwritten existing embeddings, if any?')
                case 'Yes'
                    overwrittenold = true;
                case 'No'
                    overwrittenold = false;
                otherwise
                    return;
            end
            
            [K, usehvgs] = gui.i_gethvgnum(sce, FigureHandle);
            if isempty(K), return; end                
            fw = gui.myWaitbar(FigureHandle);
            try
                forced = true;
                %if contains(methoddimtag, 'tsne'), disp('tSNE perplexity = 30'); end
                for k=1:length(methodtag)
                    if ~overwrittenold && ...
                        isfield(sce.struct_cell_embeddings, methodtag{k}) && ...
                        ~isempty(sce.struct_cell_embeddings.(methodtag{k}))
                        fprintf('Embedding cells using %s - skipping...\n', ...
                            upper(methodtag{k}));
                        continue;
                    else
                        gui.myWaitbar(FigureHandle, fw, false, '', ...
                            sprintf('Embedding cells using %s', ...
                            upper(methodtag{k})), ...
                            (k-1)/length(methodtag));
                    end
                    ndim = 2 + contains(methodtag{k},'3d');
                    sce = sce.embedcells(methodtag{k}, forced, ...
                        usehvgs, ndim, K, [], false);
                end
                % disp('Following the library-size normalization and log1p-transformation, we visualized similarity among cells by projecting them into a reduced dimensional space using t-distributed stochastic neighbor embedding (t-SNE)/uniform manifold approximation and projection (UMAP).')
                
            catch ME
                gui.myWaitbar(FigureHandle, fw);
                gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
                % rethrow(ME)
                return;
            end
            if length(methodtag)>1, readytoshow = true; end
            gui.myWaitbar(FigureHandle, fw);

        else
            return;
        end
        guidata(FigureHandle, sce);
        in_RefreshAll(src, [], true, false);   % keepview, keepcolr
        if readytoshow            
            if strcmp(gui.myQuestdlg(FigureHandle, 'Show all avaiable embeddings?'), 'Yes')
                gui.sc_multiembeddingview(sce, [], FigureHandle);
            end
        end
    end

    function in_DetermineCellTypeClustersGeneral(src, ~, usedefaultdb)
        if nargin < 3, usedefaultdb = true; end
        if usedefaultdb
            organtag = "all";
            databasetag = "panglaodb";
            if ~gui.gui_showrefinfo('PanglaoDB [PMID:30951143]', FigureHandle), return; end
            %if isempty(speciestag)
                speciestag = gui.i_selectspecies(2, false, FigureHandle);
            %end
            if isempty(speciestag), return; end
        else
            [Tm, Tw] = pkg.i_markerlist2weight(sce, FigureHandle);
            if isempty(Tm) || isempty(Tw)
                return;
            end
            wvalu = Tw.Var2;
            wgene = string(Tw.Var1);
            celltypev = string(Tm.Var1);
            markergenev = string(Tm.Var2);
        end

        [manuallyselect, bestonly] = gui.i_annotemanner(FigureHandle);

        dtp = findobj(h, 'Type', 'datatip');
        delete(dtp);
        cLdisp = cL;
        if ~manuallyselect, fw = gui.myWaitbar(FigureHandle); end

        for ix = 1:max(c)
            if ~manuallyselect
                gui.myWaitbar(FigureHandle, fw, false, '', '', ix/max(c));
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

            if gui.i_isuifig(FigureHandle)
                [indx, tf] = gui.myListdlg(FigureHandle, ctxt, 'Select cell type');
            else
                [indx, tf] = listdlg('PromptString', {'Select cell type'}, ...
                    'SelectionMode', 'single', 'ListString', ctxt, 'ListSize', [220, 300]);
            end

            if tf ~= 1, return; end
                ctxt = Tct.C1_Cell_Type{indx};
            else
                ctxt = Tct.C1_Cell_Type{1};
            end


            hold(hAx,'on');

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
            hold(hAx,'off');
        end
        if ~manuallyselect, gui.myWaitbar(FigureHandle, fw); end
        sce.c_cell_type_tx = string(cL(c));
        nx = length(unique(sce.c_cell_type_tx));
        if nx > 1
            newtx = erase(sce.c_cell_type_tx, "_{"+digitsPattern+"}");
            if length(unique(newtx)) ~= nx                
                if strcmp(gui.myQuestdlg(FigureHandle, 'Merge subclusters of same cell type?'), 'Yes')
                    in_MergeSubCellTypes(src);
                end
            end
        end
        guidata(FigureHandle, sce);
    end

    function [iscelltype] = in_pickcelltypeclusterid(a)        
        switch gui.myQuestdlg(FigureHandle, a, '', {'Cluster', 'Cell Type'}, 'Cluster')
            case 'Cluster'
                iscelltype = false;
            case 'Cell Type'
                iscelltype = true;
            otherwise
                iscelltype = [];
        end
    end

    function [hasbrushed, ptsSelected] = i_checkbrusheddata
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            hasbrushed = false;
            gui.myWarndlg(FigureHandle, "No cells are selected.");
        else
            hasbrushed = true;
        end
    end

    function in_Brushed2NewCluster(~, ~)
        [hasbrushed, ptsSelected] = i_checkbrusheddata;
        if ~hasbrushed, return; end
        [iscelltype] = in_pickcelltypeclusterid('Make a new cluster or new cell type group out of brushed cells?');
        if isempty(iscelltype), return; end
        if iscelltype
            n = sum(contains(unique(sce.c_cell_type_tx), "New cell type"));
            if n > 0
                nname = sprintf('New cell type %d', n+1);
            else
                nname = 'New cell type';
            end
            
            if gui.i_isuifig(FigureHandle)                
                newctype = gui.myInputdlg({'Enter new cell type name:'}, ...
                    '', {nname}, FigureHandle);
            else
                newctype = inputdlg('Enter new cell type name:', ...
                    '', [1, 50], {nname});
            end

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
        [h] = gui.i_gscatter3(sce.s, c, methodid, [], hAx);
        title(hAx, sce.title);
        subtitle(hAx, '[genes x cells]');
        view(hAx, [ax, bx]);
        ix_labelclusters(true);
        guidata(FigureHandle, sce);
    end

    function in_Brushed2MergeClusters(~, ~)
        [hasbrushed, ptsSelected] = i_checkbrusheddata;
        if ~hasbrushed, return; end
        [iscelltype] = in_pickcelltypeclusterid('Merge brushed cells into same cluster or same cell type?');
        if isempty(iscelltype), return; end
        if iscelltype
            c_members = unique(sce.c_cell_type_tx(ptsSelected));
        else
            c_members = unique(c(ptsSelected));
        end

        if isscalar(c_members)
            gui.myWarndlg(FigureHandle, ...
            "All brushed cells are in one cluster or belong to the same cell type.",'');
            return;
        end

        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, c_members, 'Select target cluster');
        else
            [indx, tf] = listdlg('PromptString', ...
                {'Select target cluster'}, 'SelectionMode', ...
                'single', 'ListString', string(c_members), ...
                'ListSize', [220, 300]);
        end
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
        [h] = gui.i_gscatter3(sce.s, c, methodid, [], hAx);
        title(hAx, sce.title);
        subtitle(hAx, '[genes x cells]');
        view(hAx, [ax, bx]);
        ix_labelclusters(true);
        guidata(FigureHandle, sce);
    end


    function in_Brush4Celltypes(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            gui.myHelpdlg(FigureHandle,...
                "No cells are selected. Please use the data brush tool to " + ...
                "select cells for cell type assignment.");
            return;
        end        
        if ~strcmp(gui.myQuestdlg(FigureHandle, ...
                ['This is a one-time analysis. Cell type labels will ' ...
                'not be saved. Continue?']), 'Yes')
            return;
        end
        if ~exist('speciestag', 'var') || isempty(speciestag)
            % xxx
            speciestag = gui.i_selectspecies(2, false, FigureHandle);
        end
        if isempty(speciestag), return; end
        fw = gui.myWaitbar(FigureHandle);
        [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, sce.s, ...
            ptsSelected, ...
            speciestag, "all", "panglaodb", false);
        ctxt = Tct.C1_Cell_Type;
        gui.myWaitbar(FigureHandle,fw);

        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, ctxt, 'Select cell type');
        else        
            [indx, tf] = listdlg('PromptString', ...
                {'Select cell type'}, 'SelectionMode', 'single', ...
                'ListString', ctxt, 'ListSize', [220, 300]);
        end
        if tf == 1
            ctxt = Tct.C1_Cell_Type{indx};
        else
            return;
        end
        ctxt = strrep(ctxt, '_', '\_');
        delete(findall(FigureHandle, 'Type', 'hggroup'));
        if ~exist('tmpcelltypev', 'var') || length(tmpcelltypev) ~= sce.NumCells
            tmpcelltypev = cell(sce.NumCells, 1);            
        end
        siv = sce.s(ptsSelected, :);
        si = mean(sce.s(ptsSelected, :));
        [k] = dsearchn(siv, si);
        idx = find(ptsSelected);
        tmpcelltypev{idx(k)} = ctxt;
        %assignin("base","tmpcelltypev",tmpcelltypev)
        rowx = dataTipTextRow('', tmpcelltypev);
        h.DataTipTemplate.DataTipRows = rowx;
        datatip(h, 'DataIndex', idx(k));
    end

    function in_ShowCellStates(src, ~)        
        [thisc, clabel, ~, newpickclabel] = gui.i_select1state(sce, ...
            false, false, true, false, FigureHandle);
        if isempty(thisc), return; end
        if strcmp(clabel, 'Cell Cycle Phase') && ~all(strcmp(unique(thisc), "undetermined"))
            sce.c_cell_cycle_tx = thisc;
        end        
        if strcmp(clabel, 'Workspace Variable...')
            clabel = gui.i_renamec(clabel, sce, newpickclabel, FigureHandle);
            sce.list_cell_attributes = [sce.list_cell_attributes, clabel];
            sce.list_cell_attributes = [sce.list_cell_attributes, thisc];
        end

        switch gui.i_selvariabletype(thisc, FigureHandle)
            case 'Categorical/Discrete'
                [c, cL] = grp2idx(thisc);
                sce.c = c;
                in_RefreshAll(src, [], true, false);
                target{1} = hAx;
                n = max(c);
                f = 0.5 * (n - 1) ./ n;
                f = 1 + f .* (1:2:2 * n);
                cb = colorbar(target{1}, 'Ticks', f, 'TickLabels', ...
                    strrep(cellstr(cL), '_', '\_'));
                cb.Label.String = strrep(clabel, '_', '\_');
                cbtogg = findall(FigureHandle, 'Tag', 'figToglColorbar');
                set(cbtogg,'State','on');
            case {'Numerical/Continuous','Unknown'}
                sce.c = thisc;
                switch gui.myQuestdlg(FigureHandle, 'Show scores in new figure?','', ...
                        {'Yes, new figure', 'No, current figure', 'Cancel'}, 'Yes, new figure')
                    case 'No, current figure'
                        in_RefreshAll(src, [], true, false);
                        target{1} = hAx;
                        target{2} = h;
                        if isnumeric(thisc)
                            set(target{2}, 'CData', thisc);
                        else
                            [c, cL] = grp2idx(thisc);
                            set(target{2}, 'CData', c);
                        end
                        cb = colorbar(target{1});
                        cb.Label.String = strrep(clabel, '_', '\_');
                        cbtogg = findall(FigureHandle, 'Tag', 'figToglColorbar');
                        set(cbtogg,'State','on');
                    case 'Yes, new figure'
                        switch gui.myQuestdlg(FigureHandle, 'Plot scatter plot type:','', ...
                                {'Heat','Stem'},'Heat')
                            case 'Heat'
                                figfun = @gui.i_heatscatterfig;
                            case 'Stem'
                                figfun = @gui.i_stemscatterfig;
                            otherwise
                                return;
                        end
                        figfun(sce, thisc, [], clabel, FigureHandle);
                        return;
                    otherwise
                        return;
                end
            otherwise
                return;
        end
        guidata(FigureHandle, sce);
    end

    function in_EnrichrHVGs(src, events)
        if ~gui.gui_showrefinfo('HVG Functional Analysis [PMID:31861624]', FigureHandle), return; end
        ptsSelected = logical(h.BrushData.');
        if any(ptsSelected)
            [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, ...
                sce, FigureHandle);
            if ~letdoit, return; end
            if sum(ptsSelected) < 200                
                if ~strcmp(gui.myQuestdlg(FigureHandle, sprintf('Too few cells (n = %d) selected, continue?', sum(ptsSelected))), 'Yes'), return; end
            end
            scetmp = sce.removecells(~ptsSelected);
            scetmp = scetmp.qcfilter(1000, 0.15, 15);
            gui.callback_EnrichrHVGs(src, events, scetmp);
        else            
            if ~strcmp(gui.myQuestdlg(FigureHandle, sprintf('All cells (n = %d) included, continue?', sce.NumCells),''), 'Yes'), return; end
            gui.callback_EnrichrHVGs(src, events);
        end
    end

    function in_DeleteSelectedCells(src, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            gui.myWarndlg(FigureHandle, "No cells are selected.");
            return;
        end
        [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, ...
            sce, FigureHandle);
        if ~letdoit, return; end

        answer = gui.myQuestdlg(FigureHandle, 'Delete selected or unselected cells?', '', ...
            {'Selected', 'Unselected'}, 'Selected');
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

    function in_DeleteBrushedOrUnbrushedCells(src, ~, action)
        if nargin<3, action='brushed'; end

        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            gui.myWarndlg(FigureHandle, "No cells are selected.");
            return;
        end
        [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, ...
            sce, FigureHandle);
        if ~letdoit, return; end
        
        switch lower(action)
            case 'unbrushed'
                if strcmp(gui.myQuestdlg(FigureHandle, 'Unbrushed/Unselected cells will be deleted. Continue?'),'Yes')
                    in_deletecells(src, ~ptsSelected);
                else
                    return;
                end
            case 'brushed'
                if strcmp(gui.myQuestdlg(FigureHandle, 'Brushed/Selected cells will be deleted. Continue?'),'Yes')
                    in_deletecells(src, ptsSelected);
                else
                    return;
                end
            otherwise
                return;
        end
        guidata(FigureHandle, sce);
    end

    function in_deletecells(src, ptsSelected)
        needprogressbar = false;
        if sce.NumCells > 8000, needprogressbar = true; end
        if needprogressbar, fw = gui.myWaitbar(FigureHandle); end
        sce = sce.removecells(ptsSelected);
        if needprogressbar
            gui.myWaitbar(FigureHandle, fw);
        end
        [c, cL] = grp2idx(sce.c);
        in_RefreshAll(src, [], true, true);
    end

    function in_scTenifoldNet(src, events, methodtag)
        if numel(unique(sce.c_cell_type_tx)) > 1
            if ~strcmp(gui.myQuestdlg(FigureHandle, 'This analysis is cell type-specific; however, current SCE contains multiple cell types. Continue?',''),'Yes'), return; end
        end
        switch gui.myQuestdlg(FigureHandle, 'Perform 10-fold bootstrapping of cells?','',{'Yes 🐢','No 🐇','Cancel'},'No 🐇')
            case 'No 🐇'
                if methodtag==1
                    gui.callback_scPCNet1(src, events);
                elseif methodtag==2
                    gui.callback_scTenifoldNet2lite(src, events);
                end
            case 'Yes 🐢'
                if methodtag==1
                    gui.callback_scTenifoldNet1(src, events);
                elseif methodtag==2
                    gui.callback_scTenifoldNet2(src, events);
                end
        end
    end

    function in_DrawKNNNetwork(~, ~)
        k = gui.i_inputnumk(4);
        if isempty(k), return; end
        fw = gui.myWaitbar(FigureHandle);
        set(0, 'CurrentFigure', FigureHandle);
        figure('WindowStyle', 'modal');
        sc_knngraph(sce.s, k, true);
        gui.myWaitbar(FigureHandle, fw);
    end

    function in_DrawTrajectory(src, ~)
        justload = false;
        % gui.myWarndlg(FigureHandle, 'This function should not be applied to tSNE and UMAP embeddings, as they "encourage a representation of the data as disjoint clusters, which is less meaningful for modeling continuous developmental trajectories" [PMID:25664528].', '');
        if ~isempty(f_traj) && isvalid(f_traj) && isgraphics(f_traj, 'line')
            if strcmp({f_traj.Visible}, 'on')
                switch gui.myQuestdlg(FigureHandle, 'Remove existing trajectory curve?','')
                    case 'Yes'
                        in_RefreshAll(src, [], true, true);  % keepview, keepcolr
                    case 'No'
                    otherwise
                        return;
                end
            end
        end
        if license('test', 'curve_fitting_toolbox') && ~isempty(which('cscvn'))
            answer = gui.myQuestdlg(FigureHandle, 'Which method?', '', {'splinefit', 'princurve', ...
                'manual'}, 'splinefit');
        else
            answer = gui.myQuestdlg(FigureHandle, 'Which method?', '', {'splinefit', 'princurve'}, ...
                'splinefit');
        end
        switch answer
            case 'splinefit'
                dim = 1;
                [t, xyz1] = pkg.i_pseudotime_by_splinefit(sce.s, dim, false);
                pseudotimemethod = 'splinefit';
            case 'princurve'
                [t, xyz1] = pkg.i_pseudotime_by_princurve(sce.s, false);
                pseudotimemethod = 'princurve';
            case 'manual'
                if license('test', 'curve_fitting_toolbox') && ~isempty(which('cscvn'))
                    if ~isempty(h.ZData)
                        switch gui.myQuestdlg(FigureHandle, 'This function does not work for 3D embedding. Continue to switch to 2D?')
                            case 'Yes'
                                in_Switch2D3D(src,[]);
                            otherwise
                                return;
                        end
                    end
                    if ~isempty(h.ZData), return; end

                    answer2 = gui.myQuestdlg(FigureHandle, ...
                        'Draw trajectory curve or load saved curve and pseudotime?', ...
                        '',{'Draw Curve', 'Load Saved', 'Cancel'}, 'Draw Curve');

                    switch answer2
                        case 'Load Saved'
                            [file, path] = uigetfile('*.mat', ...
                                'Select a MAT-file to Load');
                            if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end
                            if isequal(file, 0)
                                disp('User canceled the file selection.');
                                return;
                            end
                            
                            fullFileName = fullfile(path, file);
                            loadedData = load(fullFileName);
                            if isfield(loadedData, 't') && isfield(loadedData, 'xyz1') && isfield(loadedData, 'pseudotimemethod')
                                t=loadedData.t;
                                xyz1=loadedData.xyz1;
                                pseudotimemethod=loadedData.pseudotimemethod;
                            else
                                gui.myErrordlg(FigureHandle, 'Not a valid .mat file.','');
                                return;
                            end
                            justload = true;
                        case 'Draw Curve'
                            % Collect points interactively
                            % [x, y] = ginput; % 
                            if ~strcmp(gui.myQuestdlg(FigureHandle, ...
                                    ['Click on the figure to select points,' ...
                                    ' press Enter to finish.']),'Yes')
                                return; 
                            end
                            x = []; y = [];

%{
% Create uifigure and axes
fig = uifigure;
ax = uiaxes(fig);
plot(ax, rand(10,1), rand(10,1), 'o'); % Example plot
% Set a callback for mouse clicks
fig.WindowButtonDownFcn = @(src, event) disp(event.IntersectionPoint);
%}

                            hold(hAx, 'on');
                            while true
                                [xi, yi, button] = ginput(1);                        
                                % Exit the loop if Enter (ASCII 13) is pressed
                                if isempty(button) || button == 13
                                    break;
                                end
                                x = [x; xi];
                                y = [y; yi];                        
                                plot(hAx, xi, yi, 'ro', 'MarkerSize', 8, ...
                                    'LineWidth', 2);
                            end
                            hold(hAx, 'off');



                            % Fit and plot the spline curve
                            splineCurve = cscvn([x'; y']);
                            [xyz1] = fnplt(splineCurve);  % Plot the spline curve in red
                            xyz1 = xyz1';
                            
                            [t] = dsearchn(xyz1, sce.s);
                            t = (t + randn(size(t)))';
                            t = normalize(t, 'range');
                            t = t(:);
                            pseudotimemethod = 'manual';                            
                        otherwise
                            return;
                    end
                end
            otherwise
                return;
        end
        hold(hAx, 'on');        

        if size(xyz1, 2) >= 3            
            % assignin('base','xyz1',xyz1)
            if size(xyz1, 1) < 0.8*(sce.NumCells)
                xyz1 = pkg.i_interp3d(xyz1, sce.NumCells);
            end

            f_traj = plot3(hAx, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-r', 'linewidth', 2);
            t1 = text(hAx, xyz1(1, 1), xyz1(1, 2), xyz1(1, 3), 'Start', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            t2 = text(hAx, xyz1(end, 1), xyz1(end, 2), xyz1(end, 3), 'End', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        elseif size(xyz1, 2) == 2
            f_traj = plot(hAx, xyz1(:, 1), xyz1(:, 2), '-r', 'linewidth', 2);
            t1 = text(hAx, xyz1(1, 1), xyz1(1, 2), 'Start', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            t2 = text(hAx, xyz1(end, 1), xyz1(end, 2), 'End', ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        end
        hold(hAx, 'off');
        drawnow;
        pause(2);

        % pseudotimemethod
        if ~strcmp(answer, 'manual')
            switch gui.myQuestdlg(FigureHandle, 'Swap ''Start'' and ''End''?','')
                case 'Yes'
                    t1.String = 'End';
                    t2.String = 'Start';
                    t = 1 - t;
                case 'Cancel'
                    return;
            end
        end

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
        
        % -----------
        if strcmp(answer, 'manual') && ~justload
            switch gui.myQuestdlg(FigureHandle, 'Save manual trajectory curve and pseudotime to an .mat file?','')
                case 'Yes'
                    [file, path] = uiputfile('*.mat', 'Save as', ...
                        'pseudotime_manual_trajectory.mat');
                    if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end
                    if isequal(file, 0) || isequal(path, 0)
                        disp('User canceled the file selection.');
                        return;
                    end                                                                        
                    fullFileName = fullfile(path, file);
                    save(fullFileName, 'xyz1', 't', 'pseudotimemethod');
                    disp(['Variables saved to ', fullFileName]);
                case 'No'
                case 'Cancel'
                    return;
                otherwise
                    return;
            end
        end
        % --------------
        
        switch gui.myQuestdlg(FigureHandle, 'View expression of selected genes','')
            case 'Yes'
                gui.sc_pseudotimegenes(sce, t, FigureHandle);
        end
    end

    function in_ClusterCellsS(src, ~)        
        if ~strcmp(gui.myQuestdlg(FigureHandle, 'Cluster cells using embedding S?',''), 'Yes'), return; end
        [sx] = gui.i_pickembedvalues(sce,[],FigureHandle);
        if isempty(sx), return; end

        answer = gui.myQuestdlg(FigureHandle, 'Which method?', 'Select Algorithm', ...
            {'K-means 🐇', 'SnnDpc [DOI:10.1016/j.ins.2018.03.031] 🐢'}, 'K-means 🐇');
        if strcmpi(answer, 'K-means 🐇')
            methodtag = "kmeans";
        elseif strcmpi(answer, 'SnnDpc [DOI:10.1016/j.ins.2018.03.031] 🐢')
            methodtag = "snndpc";
            if ~gui.gui_showrefinfo('SnnDpc [DOI:10.1016/j.ins.2018.03.031]', FigureHandle), return; end
        else
            return;
        end
        in_reclustercells(src, methodtag, sx);
        guidata(FigureHandle, sce);
    end

    function in_ClusterCellsX(src, ~)
        if ~strcmp(gui.myQuestdlg(FigureHandle, 'Cluster cells using expression matrix X?',''), 'Yes'), return; end
        % methodtagvx = {'specter (31 secs) 🐇', 'sc3 (77 secs) 🐇', ...
        %     'simlr (400 secs) 🐢', ...
        %     'soptsc (1,182 secs) 🐢🐢', 'sinnlrr (8,307 secs) 🐢🐢🐢',};
        % methodtagv = {'specter', 'sc3', 'simlr', 'soptsc', 'sinnlrr'};
        methodtagvx = {'SC3 [PMID:28346451] 🐢🐢'};
        methodtagv = {'sc3'};
        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, methodtagvx, 'Select a clustering algorithm');
        else        
            [indx, tf] = listdlg('PromptString', ...
                {'Select a clustering algorithm'}, ...
                'SelectionMode', 'single', ...
                'ListString', methodtagvx, 'ListSize', [220, 300]);
        end
        if tf == 1
            methodtag = methodtagv{indx};
        else
            return;
        end
        if (ismcc || isdeployed)
            if strcmp(methodtag, 'sc3')
                gui.myWarndlg(FigureHandle, 'SC3 is not working in standalone application.');
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
            answer1 = gui.myQuestdlg(FigureHandle, sprintf('Using existing %s clustering?', ...
                upper(methodtag)), '', ...
                {'Yes, use existing', 'No, re-compute', ...
                'Cancel'}, 'Yes, use existing');
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
            fw = gui.myWaitbar(FigureHandle);
            try
                % [sce.c_cluster_id]=sc_cluster_x(sce.X,k,'type',methodtag);
                sce = sce.clustercells(k, methodtag, true, sx);
            catch ME
                gui.myWaitbar(FigureHandle, fw, true);
                gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
                return
            end
            gui.myWaitbar(FigureHandle, fw);
        end
        [c, cL] = grp2idx(sce.c_cluster_id);
        sce.c = c;
        in_RefreshAll(src, [], true, false);
        guidata(FigureHandle, sce);
    end

    function in_highlightcellgroups(src, ~)        
        if ~strcmp(gui.myQuestdlg(FigureHandle, 'Before highlighting cell groups, the grouping colormap will be reset. Continue?',''), 'Yes'), return; end
        in_RefreshAll(src, [], true, true);        
        if ~strcmp(gui.myQuestdlg(FigureHandle, 'Select one or more cell groups to be highlighted. Continue?',''), 'Yes'), return; end
        [idx] = gui.i_selmultidlg(cL, [], FigureHandle);
        if isempty(idx), return; end
        if idx == 0, return; end
        cm = colormap(hAx);
        cm_new = ones(size(cm))*0.9;
        if numel(cL) == height(cm)
            cm_new(idx, :) = cm(idx, :);
        end
        colormap(hAx, cm_new);
    end

    function in_labelcellgroups(src, ~)
        %a = findall(FigureHandle, 'tag', 'figMenuCellGroups___');
        %b = findall(FigureHandle, 'tag', 'figToglLabelCellGroups');
        switch src.Type            
            case 'uitoggletool'
                statetag = 'State';
            case 'uimenu'
                statetag = 'Checked';
            case 'uipushtool'
                statetag = '';
        end
        % state = src.(statetag);
        dtp = findobj(h, 'Type', 'datatip');
        % disp('...state...')
        if ~isempty(dtp) % switch from on to off
            % dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);

            if ~isempty(statetag), set(src, statetag, 'off'); end
         %   if ~isempty(a), set(a, 'Checked', 'off'); end
         %   if ~isempty(b), set(b, 'State', 'off'); end
        else
            [thisc, clabel] = gui.i_select1class(sce, true, ...
                [], [], FigureHandle);
            if isempty(thisc)
                if ~isempty(statetag), set(src, statetag, 'off'); end
                return;
            end
            [c, cL] = grp2idx(thisc);
            sce.c = c;
            in_RefreshAll(src, [], true, false);
            fprintf('Cells are colored by %s.\n', lower(clabel));
            if max(c) <= 200
                if ix_labelclusters(true)
                    if ~isempty(statetag), set(src, statetag, 'on'); end
                    %if ~isempty(a), set(a,'Checked','on'); end
                    %if ~isempty(b), set(b,'State','on'); end
                else
                    if ~isempty(statetag), set(src, statetag, 'off'); end
                    %if ~isempty(a), set(a,'Checked','off'); end
                    %if ~isempty(b), set(b,'State','off'); end
                end
            else
                if ~isempty(statetag), set(src, statetag, 'off'); end
                %if ~isempty(a), set(a,'Checked','off'); end
                %if ~isempty(b), set(b,'State','off'); end
                gui.myWarndlg(FigureHandle, 'Labels are not showing. Too many categories (n>200).');
            end
            % setappdata(FigureHandle, 'cL', cL);
            guidata(FigureHandle, sce);
            % colormap(lines(min([256 numel(unique(sce.c))])));
        end
    end

    % function [txt] = i_myupdatefcnx(pdt, ~)
    %     % pos = event_obj.Position;
    %     % idx = event_obj.DataIndex;
    %     % txt = cL(c(idx));
    %     % https://www.mathworks.com/matlabcentral/answers/549567-disable-datatips-on-click
    %     pdt.Visible = 'off';
    %     txt = '';
    % end

    function [isdone] = ix_labelclusters(notasking)
        if nargin < 1, notasking = true; end
        isdone = false;
        if ~isempty(cL)
            if notasking
                stxtyes = cL(c);
            else
                [~, cLx] = grp2idx(c);
                if isequal(cL, cLx)
                    stxtyes = c;
                else
                    stxtyes = cL(c);
                end
            end
            dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);
            if isstring(stxtyes) || iscellstr(stxtyes)
                stxtyes = strrep(stxtyes, "_", "\_");
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
                % [~, k] = medoid(siv);  geometric_median
                datatip(h, 'DataIndex', idx(kb));
            end
            % ptlabelclusters.State = 'on';
            isdone = true;
        end
    end
end


