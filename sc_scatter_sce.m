function varargout = sc_scatter_sce(sce, varargin)

if usejava('jvm') && ~feature('ShowFigureWindows')
    error('MATLAB is in a text mode. This function requires a GUI-mode.');
end
if nargin < 1
    % error('Usage: sc_scatter_sce(sce)');
    scgeatool;
    return;
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
ax = [];
bx = [];
tmpcelltypev=cell(sce.NumCells,1);

if isempty(c_in)
    sce.c = ones(size(sce.X, 2), 1);
else
    sce.c = c_in;
end
if ~isempty(s_in), sce.s = s_in; end

[c, cL] = grp2idx(sce.c);

if ~(ismcc || isdeployed)
    tagx='on';
else
    tagx='off';
end

FigureHandle = figure('Name', 'SCGEATOOL :: Single-Cell Gene Expression Analysis Tool', ...
    'position', round(1.25 * [0 0 560 420]), ...
    'visible', 'off', 'NumberTitle',tagx);
movegui(FigureHandle, 'center');

% b = uipanel(FigureHandle,'Title','B','BackgroundColor','cyan');
% b.Position = [0.18 0.40 0.30 0.35];

set(findall(FigureHandle,'ToolTipString','Link/Unlink Plot'),'Visible','Off')
set(findall(FigureHandle,'ToolTipString','Edit Plot'),'Visible','Off')
set(findall(FigureHandle,'ToolTipString','Open Property Inspector'),'Visible','Off')

%a=findall(FigureHandle,'ToolTipString','New Figure');
%a.ClickedCallback = @__;
%a=findall(f,'tag','figMenuFile');

% https://undocumentedmatlab.com/articles/customizing-standard-figure-toolbar-menubar
a=findall(FigureHandle,'tag','figMenuWindow'); delete(a);
a=findall(FigureHandle,'tag','figMenuDesktop'); delete(a);
a=findall(FigureHandle,'tag','figMenuUpdateFileNew'); delete(a);
a=findall(FigureHandle,'tag','figMenuOpen'); 
a.MenuSelectedFcn='scgeatool';
a=findall(FigureHandle,'tag','figMenuFileSaveAs'); delete(a);
a=findall(FigureHandle,'tag','figMenuFileSave'); 
a.MenuSelectedFcn=@callback_SaveX;
a=findall(FigureHandle,'tag','figMenuGenerateCode'); delete(a);

a=findall(FigureHandle,'tag','figMenuFileImportData');
a.Text='Import Data Using GEO Accession...';
a.MenuSelectedFcn=@GEOAccessionToSCE;

hAx = axes('Parent', FigureHandle);

[h] = gui.i_gscatter3(sce.s, c, methodid,1,hAx);
title(hAx,sce.title);
subtitle('[genes x cells]');

dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

mfolder = fileparts(mfilename('fullpath'));

DeftToolbarHandle = findall(FigureHandle, 'Tag','FigureToolBar');  % get the figure's toolbar handle
MainToolbarHandle = uitoolbar('Parent', FigureHandle);
set(MainToolbarHandle, 'Tag', 'MainToolBar', 'HandleVisibility', 'off', 'Visible', 'on');
UserToolbarHandle = uitoolbar('Parent', FigureHandle);
set(UserToolbarHandle, 'Tag', 'UserToolBar', 'HandleVisibility', 'off', 'Visible', 'on');

% i_addbutton(1,1,@callback_Brush4Markers,"icon-mat-filter-1-20.gif","Marker genes of brushed cells");

%i_addbutton_toggle(1,0,{@togglebtfun,@turnoffuserguiding,"icon-mat-blur-off-10.gif", ...
%    "icon-mat-blur-on-10.gif",false},"Turn on/off user onboarding toolbar");

i_addbutton_toggle(1,0,{@togglebtfun,@turnoffuserguiding, ...
    "icon-mat-unfold-more-10.gif", ...
    "icon-mat-unfold-less-10.gif",false, ...
    "Turn on/off user onboarding toolbar"});
%i_addbutton_push(1,0,@call_scgeatool,"IMG00107.GIF"," ");

i_addbutton_push(1,0,@callback_ShowGeneExpr,"list.gif","Select genes to show expression")
i_addbutton_push(1,0,@ShowCellStates,"list2.gif","Show cell state")
i_addbutton_push(1,0,@SelectCellsByQC,"plotpicker-effects.gif","Filter genes and cells")

i_addbutton_toggle(1,1,{@togglebtfun,@LabelClusters, ...
    "icon-fa-tag-10b.gif","icon-fa-tags-10b.gif", ...
    false,"Label cell groups"});
i_addbutton_push(1,0,@Brushed2NewCluster,"plotpicker-glyplot-face.gif","Add brushed cells to a new group")
i_addbutton_push(1,0,@Brushed2MergeClusters,"plotpicker-pzmap.gif","Merge brushed cells to same group")
i_addbutton_push(1,0,@RenameCellTypeBatchID,"plotpicker-scatterhist.gif","Rename cell type or batch ID");
i_addbutton_push(1,0,@call_scgeatool,"IMG00107.GIF"," ");
i_addbutton_push(1,1,@ClusterCellsS,"plotpicker-dendrogram.gif","Clustering using embedding S")
i_addbutton_push(1,0,@ClusterCellsX,"icon-mw-cluster-10.gif","Clustering using expression matrix X")
i_addbutton_push(1,1,{@DetermineCellTypeClustersGeneral,true},"plotpicker-contour.gif","Assign cell types to groups")
i_addbutton_push(1,0,@Brush4Celltypes,"brush.gif","Assign cell type to selected cells");
%i_addbutton(1,0,@ShowCellStemScatter,"IMG00067.GIF","Stem scatter plot");
i_addbutton_push(1,1,@callback_Brush4Markers,"plotpicker-kagi.gif","Marker genes of brushed cells");
i_addbutton_push(1,0,@callback_FindAllMarkers,"plotpicker-plotmatrix.gif","Marker gene heatmap");
i_addbutton_push(1,0,@call_scgeatool,"IMG00107.GIF"," ");
i_addbutton_push(1,1,@gui.callback_ShowClustersPop,"plotpicker-geoscatter.gif","Show cell clusters/groups individually");
i_addbutton_push(1,0,@gui.callback_SelectCellsByClass,"plotpicker-pointfig.gif","Select cells by class");
i_addbutton_push(1,0,@DeleteSelectedCells,"plotpicker-qqplot.gif","Delete selected cells");
i_addbutton_push(1,0,@callback_SaveX,"export.gif","Export & save data");
i_addbutton_push(1,1,@EmbeddingAgain,"plotpicker-geobubble.gif","Embedding (tSNE, UMP, PHATE)");
i_addbutton_push(1,0,@Switch2D3D,"plotpicker-image.gif","Switch 2D/3D");
i_addbutton_push(1,1,@callback_CloseAllOthers,"icon-fa-cut-10.gif","Close all other figures");
i_addbutton_push(1,0,@callback_PickPlotMarker,"plotpicker-rose.gif","Switch scatter plot marker type");
i_addbutton_push(1,0,@gui.callback_PickColorMap,"plotpicker-compass.gif","Pick new color map");
% if ~(ismcc || isdeployed)
%     i_addbutton(1,0,@callback_formatfig,"xpowerpoint.gif",'Formating Figure...');
% end
i_addbutton_push(1,0,@RefreshAll,"icon-mat-refresh-20.gif","Refresh");
i_addbutton_push(0,0,@call_scgeatool,"IMG00107.GIF"," ");
%i_addbutton(0,0,@callback_CalculateCellScores,"cellscore2.gif","Calculate cell scores from list of feature genes")
%i_addbutton(0,0,@callback_ComparePotency,"plotpicker-candle.gif","Compare differentiation potency between groups");
i_addbutton_push(0,1,@gui.callback_MultiGroupingViewer,"plotpicker-arxtimeseries.gif","Multi-grouping View...");
i_addbutton_push(0,0,@gui.callback_CrossTabulation,"plotpicker-comet.gif","Cross tabulation");
%i_addbutton_push(0,0,@call_scgeatool,"IMG00107.GIF"," ");
%i_addbutton(0,1,@callback_CompareGeneBtwCls,"plotpicker-priceandvol.gif","Compare between groups");
i_addbutton_push(0,1,@callback_CellTypeMarkerScores,"cellscore.gif","Calculate signature scores for each cell");
i_addbutton_push(0,0,@callback_CompareGeneBtwCls,"cellscore2.gif","Compare between groups");
%i_addbutton_push(0,0,@call_scgeatool,"IMG00107.GIF"," ");
i_addbutton_push(0,0,@callback_DEGene2Groups,"plotpicker-boxplot.gif","Compare 2 groups (DE analysis)");
i_addbutton_push(0,0,@callback_EnrichrHVGs,"plotpicker-andrewsplot.gif","Functional enrichment analysis with HVGs");
i_addbutton_push(0,1,@callback_BuildGeneNetwork,"noun_Network_691907.gif","Build gene regulatory network");
i_addbutton_push(0,0,@callback_CompareGeneNetwork,"noun_Deep_Learning_2424485.gif","Compare two scGRNs");
i_addbutton_push(0,1,{@gui.i_savemainfig,3},"powerpoint.gif",'Save Figure to PowerPoint File...');

gui.add_3dcamera(DeftToolbarHandle, 'AllCells');

i_addbutton_push(2,0,@turnonuserguiding,"icon-fa-thumb-tack-10.gif","Turn on user guiding toolbar");
i_addbutton_toggle(2,0,{@togglebtfun, @SelectCellsByQC, ...
    "icon-mat-filter-1-10.gif","plotpicker-effects.gif", ...
    true, "Filter genes and cells"});
i_addbutton_toggle(2,0,{@togglebtfun, @EmbeddingAgain, ...
    "icon-mat-filter-2-10.gif","plotpicker-geobubble.gif", ...
    true, "Embedding (tSNE, UMP, PHATE)"});
i_addbutton_toggle(2,0,{@togglebtfun, @ClusterCellsS, ...
    "icon-mat-filter-3-10.gif", "plotpicker-dendrogram.gif", ...
    true, "Clustering using embedding S"});
i_addbutton_toggle(2,0,{@togglebtfun, @DetermineCellTypeClustersGeneral, ...
    "icon-mat-filter-4-10.gif","plotpicker-contour.gif", true, ...
    "Assign cell types to groups"});
i_addbutton_toggle(2,0,{@togglebtfun, @callback_SaveX, ...
    "icon-mat-filter-5-10.gif", "export.gif", ...
    true, "Export & save data"});

m_vie = uimenu(FigureHandle,'Text','&Multiview','Accelerator','M');
i_addmenu(m_vie,0,@gui.callback_MultiEmbeddingViewer,'Multi-embedding View...');
i_addmenu(m_vie,0,@gui.callback_MultiGroupingViewer,'Multi-grouping View...');
i_addmenu(m_vie,0,@gui.callback_CrossTabulation,'Cross Tabulation...');

m_net = uimenu(FigureHandle,'Text','&Network','Accelerator','N');
i_addmenu(m_net,0,@callback_scPCNet1,'GRN Construction - PC Regression (w/o subsampling) [PMID:33336197] 🐢...');
i_addmenu(m_net,0,@callback_scTenifoldNet1,'GRN Construction - PC Regression (w/ subsampling) [PMID:33336197] 🐢🐢 ...');
i_addmenu(m_net,1,@callback_scTenifoldNet2lite,'GRN Comparison - scTenifoldNet (w/o subsampling) [PMID:33336197] 🐢🐢 ...');
i_addmenu(m_net,0,@callback_scTenifoldNet2,'GRN Comparison - scTenifoldNet (w/ subsampling) [PMID:33336197] 🐢🐢🐢 ...');
i_addmenu(m_net,1,@callback_scTenifoldKnk1,'Virtual Gene KO - scTenifoldKnk [PMID:35510185] 🐢🐢 ...');
i_addmenu(m_net,0,@gui.callback_VirtualKOGenKI,'Virtual Gene KO - GenKI [PMID:37246643] 🐢🐢 ...');
i_addmenu(m_net,1,@callback_scTenifoldXct,'Cell-Cell Interactions (CCIs) - scTenifoldXct [PMID:36787742] 🐢🐢 ...');
i_addmenu(m_net,0,@callback_scTenifoldXct2,'Differential CCIs - scTenifoldXct [PMID:36787742] 🐢🐢🐢 ...');

m_ext = uimenu(FigureHandle,'Text','E&xternal','Accelerator','x');
i_addmenu(m_ext,0,@gui.i_setrenv,'Check R Environment');
i_addmenu(m_ext,0,@gui.i_setpyenv,'Check Python Environment');
i_addmenu(m_ext,1,@DecontX,'Detect Ambient RNA Contamination (DecontX/R) [PMID:32138770]...');
%i_addmenu(m_ext,0,@callback_SingleRCellType,'SingleR Cell Type Annotation (SingleR/R required)...');
%i_addmenu(m_ext,0,@callback_RevelioCellCycle,'Revelio Cell Cycle Analysis (Revelio/R required)...');
% i_addmenu(m_ext,0,@callback_RunSeuratSCTransform,'Run Seurat/R SCTransform (Seurat/R required)...');
i_addmenu(m_ext,0,@RunSeuratWorkflow,'Run Seurat/R Workflow (Seurat/R) [PMID:25867923]...');
i_addmenu(m_ext,0,@callback_TrajectoryAnalysis,'Pseudotime Analysis (Monocle2/R) [PMID:28825705]...');
i_addmenu(m_ext,0,@callback_RunMonocle3,'Pseudotime Analysis (Monocle3/R) [PMID:28825705]...');
i_addmenu(m_ext,1,@callback_MELDPerturbationScore,'MELD Perturbation Score (MELD/Py) [PMID:33558698]...');
i_addmenu(m_ext,0,{@SubsampleCells,2},'Geometric Sketching (geosketch/Py) [PMID:31176620]...');
i_addmenu(m_ext,0,@HarmonyPy,'Batch Integration (Harmony/Py) [PMID:31740819]...');
i_addmenu(m_ext,0,@DoubletDetection,'Detect Doublets (Scrublet/Py) [PMID:30954476]...');
i_addmenu(m_ext,1,@callback_ExploreCellularCrosstalk,'Talklr Intercellular Crosstalk [DOI:10.1101/2020.02.01.930602]...');
i_addmenu(m_ext,0,@callback_CompareGCLBtwCls,'Differential GCL Analysis [PMID:33139959]...');
i_addmenu(m_ext,0,@callback_DiffTFActivity,'Differential TF Activity Analysis...');

m_exp = uimenu(FigureHandle,'Text','Ex&perimental','Accelerator','p');
% m_exp2 = uimenu(m_exp,'Text','sc&Tenifold Suite','Accelerator','T');
% i_addmenu(m_exp2,1,@callback_scTenifoldNet1,'scTenifoldNet - GRN Construction 🐢🐢 ...');
% i_addmenu(m_exp2,0,@callback_scTenifoldNet2,'scTenifoldNet - GRN Comparison 🐢🐢🐢 ...');
% i_addmenu(m_exp2,0,@callback_scTenifoldKnk1,'scTenifoldKnk - Virtual KO of a Gene 🐢 ...');
% i_addmenu(m_exp2,0,@callback_scTenifoldXct,'scTenifoldXct - Cell-Cell Interactions 🐢 ...');

%i_addmenu(m_exp,0,@callback_ShowPseudoTimeGenes,'Show Genes with Expression Varies with Pseudotime...');
%i_addmenu(m_exp,0,@callback_DetectCellularCrosstalk,'Ligand-Receptor Mediated Intercellular Crosstalk...');
%i_addmenu(m_exp,0,@AnnotateSubTypes,'Assign Subtypes of Cells (Neurons or T Cells)...');
i_addmenu(m_exp,0,@callback_SelectCellsByMarker,'Extract Cells by Marker(+/-) Expression...');
i_addmenu(m_exp,0,@MergeSubCellTypes,'Merge Subclusters of Same Cell Type');
%i_addmenu(m_exp,0,@callback_TwoGeneCooccurrenceTest,'Two-Gene Cooccurrence Test...');
%i_addmenu(m_exp,0,@AnnotateSubGroup,'Annotate Cell Subgroups...');
i_addmenu(m_exp,1,@WorkonSelectedGenes,'Select Top n HVGs to Work on...');
i_addmenu(m_exp,0,@SubsampleCells,'Subsample 50% Cells to Work on...');

i_addmenu(m_exp,1,@DrawKNNNetwork,'Plot Cell kNN Network...');
i_addmenu(m_exp,0,@DrawTrajectory,'Plot Cell Trajectory...');
%i_addmenu(m_exp,0,@ShowCellStemScatter,"Stem Scatter Feature Plot...");
i_addmenu(m_exp,1,@gui.callback_Violinplot,'Gene Violin Plot...');
i_addmenu(m_exp,0,@gui.callback_DrawDotplot,'Gene Dot Plot...');
i_addmenu(m_exp,0,@gui.callback_GeneHeatMap,'Gene Heatmap...');

i_addmenu(m_exp,1,@gui.callback_CalculateGeneStats,   'Calculate Gene Expression Statistics...');
i_addmenu(m_exp,0,@gui.callback_CellCycleLibrarySize, 'Library Size of Cell Cycle Phases...');
i_addmenu(m_exp,0,@gui.callback_CellCycleAssignment,  'Cell Cycle Phase Assignment...');
i_addmenu(m_exp,0,@gui.callback_ShowHgBGeneExpression,'Show HgB-genes Expression...');
i_addmenu(m_exp,0,@gui.callback_ShowMtGeneExpression, 'Show Mt-genes Expression...');
i_addmenu(m_exp,0,@gui.callback_TCellExhaustionScores,'T Cell Exhaustion Score...');
                                                  
i_addmenu(m_exp,1,{@DetermineCellTypeClustersGeneral,false},'Annotate Cell Type Using Customized Markers...');
i_addmenu(m_exp,1,{@MergeCellSubtypes,1},'Import Subtype Cell Annotation from SCE in Workspace...');
i_addmenu(m_exp,0,{@MergeCellSubtypes,2},'Import Subtype Cell Annotation from SCE Data File...');

i_addmenu(m_exp,1,@GEOAccessionToSCE,'Import Data Using GEO Accession...');
i_addmenu(m_exp,0,{@MergeSCEs,1},'Merge SCEs in Workspace...');
i_addmenu(m_exp,0,{@MergeSCEs,2},'Merge SCE Data Files...');
i_addmenu(m_exp,0,{@RenameCellTypeBatchID,'Batch ID'},'Rename Batch IDs...');
i_addmenu(m_exp,0,@callback_ViewMetaData,'View Metadata...');
i_addmenu(m_exp,1,{@gui.i_savemainfig,3},'Save Figure to PowerPoint File...');
i_addmenu(m_exp,0,{@gui.i_savemainfig,2},'Save Figure as Graphic File...');
i_addmenu(m_exp,0,{@gui.i_savemainfig,1},'Save Figure as SVG File...');
i_addmenu(m_exp,1,{@(~,~) web('https://scgeatool.github.io/')},'Visit SCGEATOOL-Standalone Website...');
i_addmenu(m_exp,0,@callback_CheckUpdates,'Check for Updates...');

% handles = guihandles( FigureHandle );
% guidata( FigureHandle, handles );

kc = numel(unique(c));
colormap(pkg.i_mycolorlines(kc));

set(FigureHandle, 'visible', 'on');
guidata(FigureHandle, sce);
% setappdata(FigureHandle,'cL',cL);

set(FigureHandle,'CloseRequestFcn',@closeRequest);

if nargout > 0
    varargout{1} = FigureHandle;
end

if ~ispref('scgeatoolbox','useronboardingtoolbar')
    gui.gui_userguidingpref(true);
    setpref('scgeatoolbox','useronboardingtoolbar',true);
else
    %if getpref('scgeatoolbox','useronboardingtoolbar')
    %    gui.gui_userguidingpref(false);
    %end
end
showuseronboarding=getpref('scgeatoolbox','useronboardingtoolbar');
if ~showuseronboarding
     set(UserToolbarHandle, 'Visible', 'off');
end


    function turnonuserguiding(~,~)
        % setpref('scgeatoolbox','useronboardingtoolbar',true);
        % set(UserToolbarHandle, 'Visible', 'on');
                    gui.gui_userguidingpref(false);

    end

    function turnoffuserguiding(~,~)
        % getpref('scgeatoolbox','useronboardingtoolbar');
        if get(UserToolbarHandle, 'Visible')=="off"
            askpref=true;
        else
            askpref=false;
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

%     function i_savefig(~,~,tag)        
%         if tag==1
%             filter = {'*.svg'};
%             [filename,filepath] = uiputfile(filter);
%             if ischar(filename)
%                 saveas(FigureHandle,[filepath filename],'svg');
%             end
%         elseif tag==2
%             % axx=gca;
%             filter = {'*.jpg';'*.png';'*.tif';'*.pdf';'*.eps'};
%             [filename,filepath] = uiputfile(filter);
%             if ischar(filename)
%                 exportgraphics(FigureHandle,[filepath filename]);
%             end
%         elseif tag==3
%             gui.i_export2pptx({FigureHandle},{'SCGEATOOL'});
%         end
%     end

    function i_addmenu(menuHdl,sepTag,callbackFnc,tooltipTxt)
        if ischar(callbackFnc) || isstring(callbackFnc)
            callbackFnc=str2func(callbackFnc);
        end
        if sepTag==1
            septag='on';
        else
            septag='off';
        end
        uimenu(menuHdl,'Text',tooltipTxt,...
            'Separator',septag,...
            'Callback',callbackFnc);
    end

    function i_addbutton_push(toolbarHdl,sepTag,callbackFnc,imgFil,tooltipTxt)
        if ischar(callbackFnc) || isstring(callbackFnc)
            callbackFnc=str2func(callbackFnc);
        end
        if sepTag==1
            septag='on';
        else
            septag='off';
        end
        if toolbarHdl==0
            barhandle=DeftToolbarHandle;
        elseif toolbarHdl==1
            barhandle=MainToolbarHandle;
        elseif toolbarHdl==2
            barhandle=UserToolbarHandle;
        end
        pt = uipushtool(barhandle, 'Separator', septag);
        pt.CData = i_get_ptImage(imgFil);
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
    end


    function i_addbutton_toggle(toolbarHdl,sepTag,callbackFnc)
        imgFil=callbackFnc{3};
        tooltipTxt=callbackFnc{6};
        %if ischar(callbackFnc{1}) || isstring(callbackFnc{1})
        %    callbackFnc=str2func(callbackFnc{1});
        %end
        if toolbarHdl==0
            barhandle=DeftToolbarHandle;
        elseif toolbarHdl==1
            barhandle=MainToolbarHandle;
        elseif toolbarHdl==2
            barhandle=UserToolbarHandle;
        end
        pt =  uitoggletool (barhandle, 'Separator', sepTag);        
        pt.CData = i_get_ptImage(imgFil);
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
    end

    function togglebtfun(src,~,func,~,imgFil, ...
            actiondelay,tooltipTxt)
        if nargin<6, actiondelay=true; end
        src.CData = i_get_ptImage(imgFil);
        if actiondelay
            if src.State=="off"
                func(src);
            else
                s='To execute the function, click the button again or locate and click the same button in the toolbar above. Hover over the button to view a description of its function.';
                uiwait(helpdlg(sprintf('%s\n%s',upper(tooltipTxt),s),''));
            end
        else
            func(src);
        end
    end

    function [ptImage]=i_get_ptImage(imgFil)
        try
            [img, map] = imread(fullfile(mfolder, 'resources', imgFil));
            ptImage = ind2rgb(img, map);
        catch
            ptImage = rand(16,16,3);
        end
    end

% ------------------------
% Callback Functions
% ------------------------
    function callback_formatfig(~,~)
        title(''); subtitle('');
        a1=xlim; b1=ylim; c1=zlim;
        hold on
        if size(sce.s, 2) > 2 && ~isempty(h.ZData)
            set(hAx,'XTick',[]); set(hAx,'YTick',[]); set(hAx,'ZTick',[]);
            hAx.XAxis.Visible = 'off'; hAx.YAxis.Visible = 'off'; hAx.ZAxis.Visible = 'off';
            quiver3(hAx,a1(1),b1(1),c1(1),0,0,20,'Color','k');
            quiver3(hAx,a1(1),b1(1),c1(1),0,20,0,'Color','k');
            quiver3(hAx,a1(1),b1(1),c1(1),20,0,0,'Color','k');
            text(a1(1)+22,b1(1),c1(1),'tSNE1');
            text(a1(1),b1(1)+22,c1(1),'tSNE2');
            text(a1(1),b1(1),c1(1)+22,'tSNE3');
            xlim(a1);ylim(b1);zlim(c1);
        else
            set(hAx,'XTick',[]); set(hAx,'YTick',[]);
            hAx.XAxis.Visible = 'off'; hAx.YAxis.Visible = 'off';            
            quiver(hAx,a1(1),b1(1),20,0,'Color','k');
            quiver(hAx,a1(1),b1(1),0,20,'Color','k');                                    
            text(a1(1)+22,b1(1),'tSNE1');
            text(a1(1),b1(1)+22,'tSNE2','Rotation',90);
            xlim(a1);ylim(b1);
        end
        
        
%       set(hAx,'xcolor',[1 1 1]); set(hAx,'ycolor',[1 1 1]); set(hAx,'zcolor',[1 1 1]);
%         if rand>0.5
%             xlabel('tSNE1'); ylabel('tSNE2'); zlabel('tSNE3');
%         else
%             xlabel('UMAP1'); ylabel('UMAP2'); zlabel('UMAP3');
%         end
        grid off
        box off
%         a1=xlim; b1=ylim; c1=zlim;
%         hold on
%         arrow3([a1(1),b1(1),c1(1)],[a1(1)+20,b1(1),c1(1)]);
%         arrow3([a1(1),b1(1),c1(1)],[a1(1),b1(1)+20,c1(1)]);
%         arrow3([a1(1),b1(1),c1(1)],[a1(1),b1(1),c1(1)+20]);
        %arrow('start',[a1(1) b1(1) c1(1)],'stop',[a1(1)+20 b1(1) c1(1)]);
        %arrow('start',[a1(1) b1(1) c1(1)],'stop',[a1(1) b1(1)+20 c1(1)]);
        %arrow('start',[a1(1) b1(1) c1(1)],'stop',[a1(1) b1(1) c1(1)+20]);
    end

    function call_scgeatool(~,~)    
        % scgeatool;
        % P = get(FigureHandle,'Position');
        % k=1;
        % set(FigureHandle,'Position',[P(1)-30*k P(2)-30*k P(3) P(4)]);
    end

    function closeRequest(hObject,~)
        if ~(ismcc || isdeployed)
%             group ='scgeatool';
%             pref = 'savedatabeforeclosing';
%             title = 'Closing SCGEATOOL';
%             quest = {'Do you want to save SCE data before closing?'
%                      ''
%                      'If you do not save the data, all changes will be lost'};
%             pbtns = {'Yes','No','Cancel'};
%             [pval,tf] = uigetpref(group,pref,title,quest,pbtns);
            ButtonName = questdlg('Save SCE before closing SCGEATOOL?');
            
            switch lower(ButtonName)
                case 'yes'
                    % tf=callback_SaveX(hObject,a);
                    if ~isempty(callinghandle)
                        guidata(callinghandle,sce);
                        delete(hObject);
                        helpdlg('SCE updated.');
                    else
                        labels = {'Save SCE to variable named:'}; 
                        vars = {'sce'};
                        sce = guidata(FigureHandle);
                        values = {sce};
                        [~,tf]=export2wsdlg(labels,vars,values,...
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

    function GEOAccessionToSCE(src,~)
        answer = questdlg('Current SCE will be replaced. Continue?');
        if ~strcmp(answer, 'Yes'), return; end
    
        acc=inputdlg({'GEO accession:'},'',[1 40],{'GSM3308545'});
        % [acc]=gui.i_inputgenelist(["GSM3308545","GSM3308546","GSM3308547"]);
        if isempty(acc), return; end
        acc=acc{1};
        if strlength(acc)>4 && ~isempty(regexp(acc,'G.+','once'))
            try                
                fw=gui.gui_waitbar;
                [sce]=sc_readgeoaccession(acc);
                [c,cL]=grp2idx(sce.c);
                gui.gui_waitbar(fw);            
                guidata(FigureHandle, sce);
                RefreshAll(src, 1, false, false);
            catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
            end
        end        
    end
    
    function MergeCellSubtypes(src, ~, sourcetag)
        [requirerefresh]=gui.callback_MergeCellSubtypes(src,[],sourcetag);
        if requirerefresh
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c_cell_type_tx);
            RefreshAll(src, 1, true);
            i_labelclusters(true);
        end
    end

    function MergeSCEs(src, ~, sourcetag)
        [requirerefresh,s]=gui.callback_MergeSCEs(src,sourcetag);
        if requirerefresh && ~isempty(s)
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c_batch_id);
            RefreshAll(src, 1, true);
            helpdlg(sprintf('%s SCEs merged.', upper(s)),'');
        end
    end

%     function AnnotateSubTypes(~, ~)
%         manuallyselect=false;
%         bestonly=true;
%         speciestag = gui.i_selectspecies;
%         if isempty(speciestag), return; end
%         subtype = questdlg('Which cell type?',...
%             'Select Cell Type', 'Neurons', 'T Cells','Neurons');
%         if isempty(subtype), return; end
%         organtag = "all";
%         databasetag = "panglaodb";
%         dtp = findobj(h,'Type','datatip');
%         delete(dtp);
%         cLdisp = cL;
%         if ~manuallyselect, fw=gui.gui_waitbar_adv; end
%         for i = 1:max(c)
%             if ~manuallyselect
%                 gui.gui_waitbar_adv(fw,i/max(c));
%             end
%             ptsSelected = c == i;
%             [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, ...
%                 sce.s, ptsSelected, ...
%                 speciestag, organtag, databasetag, bestonly, subtype);
%             if isempty(Tct)
%                 ctxt={'Unknown'};
%             else
%                 ctxt=Tct.C1_Cell_Type;
%             end            
%             if ~manuallyselect, gui.gui_waitbar_adv(fw); end
%         end
%     end

%     function AnnotateSubGroup(src, ~)
%         [requirerefresh,highlightindex,newclassidenty]=gui.callback_AnnotateSubGroup(src);
%         if requirerefresh && ~isempty(newclassidenty)
%             disp('OK');
%             sce.c(highlightindex)=sce.c(highlightindex)+newclassidenty./10;
%             [c, cL] = grp2idx(sce.c);
%             RefreshAll(src, 1, true);
%             answer=questdlg('Save current classes to SCE.C_CLUSTER_ID?');
%             if strcmp(answer,'Yes')
%                 sce.c_cluster_id=c;
%                 guidata(FigureHandle,sce);
%             end
%         end
%     end

    function WorkonSelectedGenes(src,~)
%         answer=questdlg('Input the number of HVGs. Continue?');
%         if ~strcmp(answer,'Yes'), return; end
        k = gui.i_inputnumk(2000,1,sce.NumGenes,'the number of HVGs');
        if isempty(k), return; end
        answer = questdlg('Which method?','Select Method', ...        
            'Brennecke et al. (2013)','Splinefit Method',...
            'Brennecke et al. (2013)');
        fw = gui.gui_waitbar;
        switch answer
            case 'Brennecke et al. (2013)'            
                T=sc_hvg(sce.X,sce.g);
            case 'Splinefit Method'
                T=sc_splinefit(sce.X,sce.g);
            otherwise
                return;
        end
        glist=T.genes(1:min([k, sce.NumGenes]));
        % [glist]=gui.i_selectngenes(sce);
        % glist=unique(glist);
        [y,idx]=ismember(glist,sce.g);
        if ~all(y), errordlg('Runtime error.'); return; end
        sce.g=sce.g(idx);
        sce.X=sce.X(idx,:);
        gui.gui_waitbar(fw);
        RefreshAll(src, 1, true);
    end

    function GeoSketching(scr,~)
        helpdlg('Under Development.');
    end

    function SubsampleCells(src,~,methodoption)
        if nargin<3, methodoption=[]; end
        answer1=questdlg('This function subsamples 50% of cells. Continue?');
        if ~strcmp(answer1,'Yes'), return; end
        
        if isempty(methodoption)
        answer=questdlg('Select method:','', ...
            'Uniform Sampling', ...
            'Geometric Sketching [PMID:31176620]','Uniform Sampling');
        switch answer
            case 'Uniform Sampling'
                methodoption=1;
            case 'Geometric Sketching [PMID:31176620]'
                methodoption=2;
            otherwise
                return;
        end
        end

        %fw=gui.gui_waitbar;
        tn=round(sce.NumCells/2);        
        if methodoption==1
            idx=randperm(sce.NumCells);
            ids=idx(1:tn);
        elseif  methodoption==2
            gui.gui_showrefinfo('Geometric Sketching [PMID:31176620]');
            answerx=questdlg('This method requires Python environment and geosketch package installed. Continue?');
            if ~strcmp(answerx,'Yes'), return; end
            Xn=log(1+sc_norm(sce.X))';
            
            try
                [~,Xn]=pca(Xn,'NumComponents',300);
                ids=run.py_geosketch(Xn,tn);
            catch ME
                %gui.gui_waitbar(fw,true);
                errordlg(ME.message);
                return;
            end
        end
        if ~isempty(ids)
            sce=sce.selectcells(ids);
            c=sce.c;
            %gui.gui_waitbar(fw);
            RefreshAll(src, 1, true);
        else
            errordlg('Running error. No action is taken.');
        end
    end

    function SelectCellsByQC(src, ~)
        oldn=sce.NumCells;
        oldm=sce.NumGenes;
        %try
            [requirerefresh,highlightindex]=...
                gui.callback_SelectCellsByQC(src);
        %catch ME
        %    errordlg(ME.message);
        %    return;
        %end
        sce = guidata(FigureHandle);
        if requirerefresh 
            [c, cL] = grp2idx(sce.c);
            RefreshAll(src, 1, true);
            newn=sce.NumCells;
            newm=sce.NumGenes;
            helpdlg(sprintf('%d cells removed; %d genes removed.',...
                oldn-newn,oldm-newm),'');
        end
        if ~isempty(highlightindex)
            h.BrushData=highlightindex;
        end
    end

%     function RunSeuratSCTransform(src,~)
%         if callback_RunSeuratSCTransform(src)
%             guidata(FigureHandle,sce);
%         end
%     end

    function RunSeuratWorkflow(src,~)
       [ok]=gui.i_confirmscript('Run Seurat/R Workflow (Seurat)?', ...
            'R_Seurat','r');
        if ~ok, return; end
       
       [ndim]=gui.i_choose2d3d;
       if isempty(ndim), return; end
	   fw = gui.gui_waitbar;
       try
           [sce]=run.r_seurat(sce,ndim);
           [c, cL] = grp2idx(sce.c);
       catch
       	   gui.gui_waitbar(fw);
           return;
       end
	   gui.gui_waitbar(fw);
       RefreshAll(src, 1, true, false);
    end

    function DecontX(~,~)
        gui.gui_showrefinfo('DecontX [PMID:32138770]');
                    
        [ok]=gui.i_confirmscript('Detect Ambient RNA Contamination (decontX)', ...
            'R_decontX','r');
        if ~ok, return; end
        fw = gui.gui_waitbar;
        try
            [Xdecon,contamination]=run.decontX(sce);
        catch
            gui.gui_waitbar(fw);
            errordlg('Runtime error.')
            return;
        end
        gui.gui_waitbar(fw);
        figure('WindowStyle','modal');
        gui.i_stemscatter(sce.s,contamination);
        % zlim([0 1]);
        zlabel('Contamination rate')
        title('Ambient RNA contamination')
        answer=questdlg("Remove contamination?");
        switch answer
            case 'Yes'
                sce.X=round(Xdecon);
                guidata(FigureHandle,sce);
                helpdlg('Contamination removed.','');
       end
    end

    function HarmonyPy(src, ~)
        gui.gui_showrefinfo('Harmony [PMID:31740819]');        
        if numel(unique(sce.c_batch_id))<2
            warndlg('No batch effect (SCE.C_BATCH_ID is empty)');
            return;
        end        
        [c1]=grp2idx(sce.c);
        [c2]=grp2idx(sce.c_batch_id);
        if ~isequal(c1,c2)
            answer=questdlg('Color cells by batch id (SCE.C_BATCH_ID)?','');
            switch answer
                case 'Yes'
                    [c,cL] = grp2idx(sce.c_batch_id);
                    sce.c = c;
                    RefreshAll(src, 1, true, false);
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
            RefreshAll(src, 1, true, false);

            ButtonName = questdlg('Update Saved Embedding?','');
            switch ButtonName
                case 'Yes'
                    [methodtag]=gui.i_pickembedmethod;
                    if isempty(methodtag), return; end        
                    if ismember(methodtag,{'tsne','umap','phate','metaviz'})
                        sce.struct_cell_embeddings.(methodtag)=sce.s;
                    end
                    helpdlg(sprintf('%s Embedding is updated.',methodtag),'');                    
            end
        end
        guidata(FigureHandle, sce);
    end


    function DoubletDetection(src, ~)
        gui.gui_showrefinfo('Scrublet [PMID:30954476]');
        [isDoublet,doubletscore,methodtag,done]=gui.callback_DoubletDetection(src);
        if done && ~any(isDoublet)
            helpdlg('No doublet detected.','');
            return;
        end
        if done && any(isDoublet) && sce.NumCells==length(doubletscore)
            tmpf_doubletdetection=figure('WindowStyle','modal');
            gui.i_stemscatter(sce.s,doubletscore);
            zlabel('Doublet Score')
            title(sprintf('Doublet Detection (%s)',methodtag))
            answer=questdlg(sprintf("Remove %d doublets?",sum(isDoublet)));
                switch answer
                    case 'Yes'
                        close(tmpf_doubletdetection);
                        % i_deletecells(isDoublet);
                        sce = sce.removecells(isDoublet);
                        guidata(FigureHandle,sce);
                        [c, cL] = grp2idx(sce.c);
                        RefreshAll(src, 1, true, false);
                        helpdlg('Doublets deleted.','');
                end
        end
    end        


    function MergeSubCellTypes(src,~)
        if isempty(sce.c_cell_type_tx), return; end
        % [sce]=pkg.i_mergeSubCellNames(sce);        
        newtx=erase(sce.c_cell_type_tx,"_{"+digitsPattern+"}");
        if isequal(sce.c_cell_type_tx,newtx)
            helpdlg("No sub-clusters are meraged.");
        else
            sce.c_cell_type_tx=newtx;
            [c,cL]=grp2idx(sce.c_cell_type_tx);
            sce.c = c;
            RefreshAll(src, 1, true, false);
            i_labelclusters;
        end
        guidata(FigureHandle, sce);
    end

% =========================
    function RefreshAll(src, ~, keepview, keepcolr)
        if nargin < 4, keepcolr = false; end
        if nargin < 3, keepview = false; end        
        if keepview || keepcolr
            [para] = gui.i_getoldsettings(src);
        end
        if size(sce.s, 2) > 2 && ~isempty(h.ZData)
            if keepview, [ax, bx] = view(); end
            h = gui.i_gscatter3(sce.s, c, methodid, hAx);
            if keepview, view(ax, bx); end
        else   % otherwise 2D
            h = gui.i_gscatter3(sce.s(:, 1:2), c, methodid, hAx);
        end
        if keepview
            h.Marker = para.oldMarker;
            h.SizeData = para.oldSizeData;
        end
        if keepcolr
            colormap(para.oldColorMap);
        else
            kc = numel(unique(c));
            colormap(pkg.i_mycolorlines(kc));
        end
        title(sce.title);
        subtitle('[genes x cells]');
        % ptlabelclusters.State = 'off';
        % UitoolbarHandle.Visible='off';
        % UitoolbarHandle.Visible='on';
        guidata(FigureHandle, sce);
    end

    function Switch2D3D(src, ~)
        [para] = gui.i_getoldsettings(src);
        if isempty(h.ZData)   % current 2 D
            if ~(size(sce.s, 2) > 2)
                helpdlg('Canno swith to 3-D. SCE.S is 2-D','');
                return;
            end
            h = gui.i_gscatter3(sce.s, c, methodid, hAx);
            if ~isempty(ax) && ~isempty(bx) && ~any([ax bx] == 0)
                view(ax, bx);
            else
                view(3);
            end
        else                 % current 3D do following
            [ax, bx] = view();
            answer = questdlg('Which view to be used to project cells?', '', ...
                'X-Y Plane', 'Screen/Camera', 'PCA-rotated', 'X-Y Plane');
            switch answer
                case 'X-Y Plane'
                    sx=sce.s;
                case 'Screen/Camera'
                    sx = pkg.i_3d2d(sce.s, ax, bx);                   
                case {'PCA-rotated'}
                    [~,sx]=pca(sce.s);
                otherwise
                    return;
            end
            h = gui.i_gscatter3(sx(:, 1:2), c, methodid, hAx);
            sce.s=sx;
        end
        title(sce.title);
        subtitle('[genes x cells]');

        h.Marker = para.oldMarker;
        h.SizeData = para.oldSizeData;
        colormap(para.oldColorMap);
    end

    function RenameCellTypeBatchID(src,~,answer)
        if nargin<3 || isempty(answer)
            answer=questdlg('Rename cell type, batch ID, or gene name?','','Cell type','Batch ID','Gene name','Cell type');
        end
        switch answer
            case 'Cell type'
                [requirerefresh]=gui.callback_RenameCellType(src);
            case 'Batch ID'
                [requirerefresh]=gui.callback_RenameBatchID(src);
            case 'Gene name'
                [requirerefresh]=gui.callback_RenameGenes(src);
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
           i_labelclusters(false);
        end
    end

    function EmbeddingAgain(src, ~)
        [methodtag]=gui.i_pickembedmethod;
        if isempty(methodtag), return; end
        if isempty(sce.struct_cell_embeddings)
            sce.struct_cell_embeddings = struct('tsne', [], 'umap', [], ...
                'phate', [], 'metaviz', []);
        end
        %methodtag = lower(answer);
        usingold = false;
        if ~isfield(sce.struct_cell_embeddings,'metaviz')
            sce.struct_cell_embeddings.('metaviz')=[];
        end
        
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
        whitelist=[];
        if ~usingold
            answer2 = questdlg(sprintf('Use highly variable genes (HVGs, n=2000) or use all genes (n=%d)?', sce.NumGenes), ...
                '', '2000 HVGs 🐇', 'All Genes 🐢', 'Other...', '2000 HVGs 🐇');
            switch answer2
                case 'All Genes 🐢'
                    usehvgs = false;
                    K=sce.NumGenes;
                case '2000 HVGs 🐇'                    
                    usehvgs = true;
                    K=2000;
                case 'Other...'
                     K=gui.i_inputnumk(min([3000,sce.NumGenes]), ...
                         100,sce.NumGenes);
                     if isempty(K), return; end
                     usehvgs = true;
                     [whitelist]=gui.i_selectwhitelist(sce);
                    if isnumeric(whitelist) 
                        if whitelist==0
                            return;
                        end
                    end                       
                otherwise
                    return;
            end
            [ndim]=gui.i_choose2d3d;
            if isempty(ndim), return; end
            if ~strcmpi(methodtag,'metaviz')
                fw = gui.gui_waitbar;
            end
            try
                forced = true;
                if strcmpi(methodtag,'tsne'), disp('tSNE perplexity = 30'); end
                sce = sce.embedcells(methodtag, forced, usehvgs, ndim, K, whitelist);
            catch ME
                if ~strcmpi(methodtag,'metaviz')
                    gui.gui_waitbar(fw,true);
                end
                errordlg(ME.message);
                return
            end
            if ~strcmpi(methodtag,'metaviz'), gui.gui_waitbar(fw); end
        end
        RefreshAll(src, 1, true, false);
        guidata(FigureHandle, sce);
        disp('Following the library-size normalization and log1p-transformation, we visualized similarity among cells by projecting them into a reduced dimensional space using t-distributed stochastic neighbor embedding (t-SNE)/uniform manifold approximation and projection (UMAP).')
    end

%     function [Tct]=i_determinecelltype(sce, ptsSelected, wvalu, wgene, celltypev, markergenev)
%             T=table(celltypev);
%             Xk=sce.X(:,ptsSelected);
%             S=zeros(length(celltypev),1);
%             for j=1:length(celltypev)
%                 g=strsplit(markergenev(j),',');
%                 g=g(1:end-1);
%                 Z=0; ng=0;
%                 for ix=1:length(g)
%                     if any(g(ix)==wgene) && any(g(ix)==sce.g)
%                         wi=wvalu(g(ix)==wgene);  
%                         z=median(Xk(sce.g==g(ix),:));
%                         Z=Z+z*wi;
%                         ng=ng+1;
%                     end
%                 end
%                 if Z>0, S(j)=Z./nthroot(ng,3); end
%             end
%             if all(S(:)==0)
%                 Tct=cell2table({'Unknown',0});
%             else
%                 [~,idx]=sort(S,'descend');
%                 T=[T,array2table(S)];
%                 Tct=T(idx,:);                
%             end
%             Tct.Properties.VariableNames={'C1_Cell_Type','C1_CTA_Score'};
%     end


    function DetermineCellTypeClustersGeneral(src, ~, usedefaultdb)
        if nargin<3, usedefaultdb=true; end
        if usedefaultdb
            organtag = "all";
            databasetag = "panglaodb";
            gui.gui_showrefinfo('PanglaoDB [PMID:30951143]');
            speciestag = gui.i_selectspecies(2);
            if isempty(speciestag), return; end
        else
            [Tm,Tw]=pkg.i_markerlist2weight(sce);
            if isempty(Tm)||isempty(Tw)
                return;
            end        
            wvalu=Tw.Var2;
            wgene=string(Tw.Var1);
            celltypev=string(Tm.Var1);
            markergenev=string(Tm.Var2);
        end

        [manuallyselect,bestonly]=i_annotemanner;

        dtp = findobj(h,'Type','datatip');
        delete(dtp);
        cLdisp = cL;
        if ~manuallyselect, fw=gui.gui_waitbar_adv; end

        for i = 1:max(c)
            if ~manuallyselect
                gui.gui_waitbar_adv(fw,i/max(c));
            end
            ptsSelected = c == i;

        if usedefaultdb
            [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, ...
                sce.s, ptsSelected, ...
                speciestag, organtag, databasetag, bestonly);
        else
            [Tct]=pkg.e_determinecelltype(sce, ptsSelected, wvalu, ...
                wgene, celltypev, markergenev);
        end

            ctxt=Tct.C1_Cell_Type;
            
            if manuallyselect && length(ctxt)>1
                [indx, tf] = listdlg('PromptString', {'Select cell type'},...
                    'SelectionMode', 'single', 'ListString', ctxt);
                if tf ~= 1, return; end
                ctxt = Tct.C1_Cell_Type{indx};
            else
                ctxt = Tct.C1_Cell_Type{1};
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
                [k] = dsearchn(siv, si);        % Nearest point search
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
        nx=length(unique(sce.c_cell_type_tx));
        if nx>1
            newtx=erase(sce.c_cell_type_tx,"_{"+digitsPattern+"}");
            if length(unique(newtx))~=nx
                answer = questdlg('Merge subclusters of same cell type?');
                if strcmp(answer, 'Yes')
                    MergeSubCellTypes(src);
                end
            end
        end
        guidata(FigureHandle, sce);
    end

    function [iscelltype]=pick_celltype_clusterid(a)
        answer = questdlg(a,'',...
                    'Cluster','Cell Type','Cluster');        
        switch answer
            case 'Cluster'
                iscelltype=false;
            case 'Cell Type'
                iscelltype=true;
            otherwise
                iscelltype=[];
        end
    end

    function Brushed2NewCluster(~, ~)
        [iscelltype]=pick_celltype_clusterid('Make a new cluster or new cell type group out of brushed cells?');
        if isempty(iscelltype), return; end
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return;
        end
        if iscelltype
            n=sum(contains(unique(sce.c_cell_type_tx),"New cell type"));
            if n>0
                nname=sprintf('New cell type %d',n+1);
            else
                nname='New cell type';
            end
            newctype = inputdlg('Enter new cell type name:', '', [1 50], {nname});
            if isempty(newctype), return; end       
            sce.c_cell_type_tx(ptsSelected) = string(newctype);
            [c, cL] = grp2idx(sce.c_cell_type_tx);            
        else
            c(ptsSelected) = max(c) + 1;
            [c, cL] = grp2idx(c);            
            sce.c_cluster_id = c;
        end
            sce.c = c;
            [ax, bx] = view();
            [h] = gui.i_gscatter3(sce.s, c, methodid, hAx);
            title(sce.title);
            subtitle('[genes x cells]');
            view(ax, bx);
            i_labelclusters(true);
            guidata(FigureHandle, sce);
    end

    function Brushed2MergeClusters(~, ~)
        [iscelltype]=pick_celltype_clusterid('Merge brushed cells into same cluster or same cell type?');
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
        
        [indx, tf] = listdlg('PromptString',...
            {'Select target cluster'}, 'SelectionMode',...
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
        [ax, bx] = view();
        [h] = gui.i_gscatter3(sce.s, c, methodid, hAx);
        title(sce.title);
        subtitle('[genes x cells]');
        view(ax, bx);
        i_labelclusters(true);
        guidata(FigureHandle, sce);
    end

    function Brush4Celltypes(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            helpdlg("No cells are selected. Please use the data brush tool to select cells for cell type assignment.",'');
            return;
        end

        answer = questdlg('This is a one-time analysis. Cell type labels will not be saved. Continue?');
        if ~strcmp(answer, 'Yes')
            return;
        end

        speciestag = gui.i_selectspecies;
        if isempty(speciestag), return; end  
        fw = gui.gui_waitbar;
        [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, sce.s,...
            ptsSelected, ...
            speciestag, "all", "panglaodb",false);
        ctxt = Tct.C1_Cell_Type;
        gui.gui_waitbar(fw);
        
        [indx, tf] = listdlg('PromptString',...
            {'Select cell type'}, 'SelectionMode', 'single', 'ListString', ctxt);
        if tf == 1
            ctxt = Tct.C1_Cell_Type{indx};
        else
            return;
        end
        ctxt = strrep(ctxt, '_', '\_');
        delete(findall(FigureHandle,'Type','hggroup'));
        if ~exist('tmpcelltypev','var')
            tmpcelltypev=cell(sce.NumCells,1);
        end
        siv=sce.s(ptsSelected, :);
        si=mean(sce.s(ptsSelected, :));
        [k] = dsearchn(siv, si);
        idx = find(ptsSelected);
        tmpcelltypev{idx(k)}=ctxt;
        row = dataTipTextRow('', tmpcelltypev);
        h.DataTipTemplate.DataTipRows = row;        
        datatip(h, 'DataIndex', idx(k));
    end        
    
    function ShowCellStemScatter(~, ~)
        sce = guidata(FigureHandle);
        [thisc,clable]=gui.i_select1state(sce);
        if isempty(thisc), return; end
        figure('WindowStyle','modal');
        gui.i_stemscatter(sce.s,grp2idx(thisc));
        zlabel(clable);        
    end

    function ShowCellStates(src, ~)
        sce=guidata(FigureHandle);        
        [thisc,clable,~,newpickclable]=gui.i_select1state(sce);        
        %clable
        %newpickclable
        if strcmp(clable,'Cell Cycle Phase')
            if length(unique(thisc))>1
                sce.c_cell_cycle_tx=thisc;
            end
        end
        if isempty(thisc), return; end
        if strcmp(clable,'Workspace Variable...')
            clable=gui.i_renamec(clable,sce,newpickclable);
            sce.list_cell_attributes=[sce.list_cell_attributes,clable];
            sce.list_cell_attributes=[sce.list_cell_attributes,thisc];
        end
        [c,cL]=grp2idx(thisc);
        sce.c=c;
        [answer]=gui.i_selvariabletype(thisc);
        RefreshAll(src, 1, true, false);
        switch answer
            case 'Categorical/Discrete'
                n=max(c);
                f=0.5*(n-1)./n;
                f=1+f.*(1:2:2*n);
                cb=colorbar('Ticks',f,'TickLabels', ...
                    strrep(cellstr(cL),'_','\_'));
            otherwise   % case 'Numerical/Continuous'
                if isnumeric(thisc)
                    set(h,'CData',thisc);
                else
                    set(h,'CData',c);
                end
                cb=colorbar;
        end
        cb.Label.String = strrep(clable,'_','\_');
        % helpdlg(clable)
        % guidata(FigureHandle, sce);
    end

    function DeleteSelectedCells(src, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.",'');
            return;
        end
        [ptsSelected,letdoit]=gui.i_expandbrushed(ptsSelected,sce);
        if ~letdoit, return; end
        
        answer = questdlg('Delete selected or unselected cells?','', ...
            'Selected', 'Unselected','Selected');
        if isempty(answer), return; end
        if strcmp(answer, 'Unselected')
            i_deletecells(src,~ptsSelected);
        elseif strcmp(answer, 'Selected')            
            i_deletecells(src,ptsSelected);
        else
            return;
        end
        guidata(FigureHandle,sce);
    end

    function i_deletecells(src,ptsSelected)
        needprogressbar=false;
        if sce.NumCells>8000, needprogressbar=true; end
        if needprogressbar
            fw = gui.gui_waitbar;
        end
        sce = sce.removecells(ptsSelected);
        if needprogressbar
            gui.gui_waitbar(fw);
        end
        [c, cL] = grp2idx(sce.c);
        RefreshAll(src,1,true,true);
    end

    function DrawKNNNetwork(~,~)
        k=gui.i_inputnumk(3);
        if isempty(k), return; end
        fw = gui.gui_waitbar;
        set(0, 'CurrentFigure', FigureHandle);
        figure('WindowStyle','modal');
        sc_knngraph(sce.s,k,true);
        gui.gui_waitbar(fw);
    end

    function DrawTrajectory(~, ~)
        warndlg('This function should not be applied to tSNE and UMAP embeddings, as they "encourage a representation of the data as disjoint clusters, which is less meaningful for modeling continuous developmental trajectories" [PMID:25664528].','');
        answer = questdlg('Which method?', '', ...
            'splinefit (🐇)', 'princurve (🐢)', ...
            'splinefit (🐇)');
        if strcmp(answer, 'splinefit (🐇)')
            dim = 1;
            [t, xyz1] = pkg.i_pseudotime_by_splinefit(sce.s, dim, false);
        elseif strcmp(answer, 'princurve (🐢)')
            [t, xyz1] = pkg.i_pseudotime_by_princurve(sce.s, false);
        else            
            return;
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
                try
                    psf1=figure('WindowStyle','modal');
                    pkg.i_plot_pseudotimeseries(log2(sce.X + 1), ...
                        sce.g, t, selectedg);
                catch ME
                    if exist('psf1','var') && ishandle(psf1)
                        close(psf1);
                    end
                    errordlg(ME.message);
                end
            case 'No'
                return;
        end        
    end

    function ClusterCellsS(src, ~)
        answer = questdlg('Cluster cells?');
        if ~strcmp(answer, 'Yes'), return; end        
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
        i_reclustercells(src, methodtag);
        guidata(FigureHandle, sce);
    end

    function ClusterCellsX(src, ~)
        answer = questdlg('Cluster cells using X?');
        if ~strcmp(answer, 'Yes'), return; end
        
        methodtagvx = {'specter (31 secs) 🐇','sc3 (77 secs) 🐇',...
             'simlr (400 secs) 🐢',...
             'soptsc (1,182 secs) 🐢🐢', 'sinnlrr (8,307 secs) 🐢🐢🐢', };
        methodtagv = {'specter','sc3','simlr', 'soptsc', 'sinnlrr'};
        [indx, tf] = listdlg('PromptString',...
            {'Select clustering program'},...
            'SelectionMode', 'single', ...
            'ListString', methodtagvx);
        if tf == 1
            methodtag = methodtagv{indx};
        else
            return;
        end
        if (ismcc || isdeployed)
            if strcmp(methodtag,'sc3')
                warndlg('SC3 is not working in standalone application.','');
                return;
            end
        end
        i_reclustercells(src, methodtag);
        guidata(FigureHandle, sce);
    end

    function i_reclustercells(src, methodtag)
        methodtag = lower(methodtag);
        usingold = false;
        if ~isempty(sce.struct_cell_clusterings.(methodtag))
            answer1 = questdlg(sprintf('Using existing %s clustering?',...
                upper(methodtag)), '', ...
                'Yes, use existing', 'No, re-compute',...
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
            k = gui.i_inputnumk;
            if isempty(k), return; end
            fw = gui.gui_waitbar;
            try
                % [sce.c_cluster_id]=sc_cluster_x(sce.X,k,'type',methodtag);
                sce = sce.clustercells(k, methodtag, true);
            catch ME
                gui.gui_waitbar(fw,true);
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
        dtp = findobj(h, 'Type', 'datatip');
        %disp('...state...')
        if strcmp(state, 'off') || ~isempty(dtp)  % switch from on to off
            %dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);
            set(src, 'State', 'off');
        else
            sce=guidata(FigureHandle);
            [thisc,clable]=gui.i_select1class(sce);
            if isempty(thisc)
                set(src, 'State', 'off');
                return; 
            end
            [c,cL] = grp2idx(thisc);
            sce.c = c;
            RefreshAll(src, 1, true, false);
            fprintf('Cells are colored by %s.\n', lower(clable));
            if max(c)<=200
                if i_labelclusters
                    set(src, 'State', 'on');
                else
                    set(src, 'State', 'off');
                end
            else
                set(src, 'State', 'off');
                warndlg('Labels are not showing. Too many categories (n>200).');
            end
            setappdata(FigureHandle,'cL',cL);
            guidata(FigureHandle, sce);              
            % colormap(lines(min([256 numel(unique(sce.c))])));
        end
    end

    function [txt] = i_myupdatefcnx(~, event_obj)
        % pos = event_obj.Position;
        idx = event_obj.DataIndex;
        txt = cL(c(idx));        
    end

    function [isdone] = i_labelclusters(notasking)
        if nargin < 1, notasking = false; end
        isdone = false;
        if ~isempty(cL)
            if notasking
                %stxtyes = c;
                stxtyes = cL(c);
            else
                [~,cLx]=grp2idx(c);
                if isequal(cL,cLx)
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
                stxtyes=strrep(stxtyes,"_","\_");
            end
            row = dataTipTextRow('', stxtyes);
            h.DataTipTemplate.DataTipRows = row;
            % h.DataTipTemplate.FontSize = 5;
            for i = 1:max(c)
                idx = find(c == i);
                siv = sce.s(idx, :);
                si = mean(siv, 1);
                % si=geometric_median(siv');
                [k] = dsearchn(siv, si);
                %[~, k] = medoid(siv);  geometric_median
                datatip(h, 'DataIndex', idx(k));
            end
            %ptlabelclusters.State = 'on';
            isdone = true;
        end
    end

end

