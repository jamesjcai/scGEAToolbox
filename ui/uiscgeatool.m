function uiscgeatool

useexample = rand>0.5;

YOffset = 5;
p = get(0,"ScreenSize");
maxh = p(4) - 2*YOffset - 50;

import pkg.*
import gui.*

mfolder = fileparts(mfilename('fullpath'));


c=[];
cL=[];
ax=[]; bx=[]; f_traj=[];

FigureHandle = uifigure('Visible', 'off');
FigureHandle.Position = round(1.25*[0, 0, 560, 420]);
FigureHandle.Name = 'scgeatool ui';
hAx = uiaxes(FigureHandle,'Visible','off');


height = 480;
width = 580;
sz = FigureHandle.Position;
%x = mean( sz( [1, 3]));
%y = mean( sz( [2, 4]));
hAx.Position= [(sz(3) - width)/2, (sz(4) - height)/2, width, height];
movegui(FigureHandle, 'center');
drawnow;

if useexample
    FigureHandle.Visible="on";
    pD = uiprogressdlg(FigureHandle,...
        Indeterminate = "on",...
        Message = "Loading...",...
        Title = "Please wait");
    load(fullfile(mfolder,'..','example_data','workshop_example.mat'),'sce');
else
    sce = SingleCellExperiment;
end

if FigureHandle.Position(4) + YOffset > maxh
    FigureHandle.Position(4) = maxh;
    scroll(FigureHandle,"top");
end

%[h]=scatter3(hAx,randn(300,1),randn(300,1),randn(300,1));
%h.Visible="off";
h=[];

%sce = SingleCellExperiment;


% dt = datacursormode;
% dt.UpdateFcn = {@i_myupdatefcnx};


DeftToolbarHandle = uitoolbar(FigureHandle);
MainToolbarHandle = uitoolbar(FigureHandle);
in_addbuttonpush(1, 0, @in_callback_ShowGeneExpr, "list.gif", "Select genes to show expression")
in_addbuttonpush(1, 0, @in_ShowCellStates, "list2.gif", "Show cell state")
in_addbuttonpush(1, 0, @in_SelectCellsByQC, "plotpicker-effects.gif", "Filter genes and cells")
in_addbuttontoggle(1, 1, {@in_togglebtfun, @in_labelcellgroups, ...
     "icon-fa-tag-10b.gif", "icon-fa-tags-10b.gif", ...
     false, "Label cell groups"});
in_addbuttonpush(1, 0, @in_Brushed2NewCluster, "plotpicker-glyplot-face.gif", "Add brushed cells to a new group")
in_addbuttonpush(1, 0, @in_Brushed2MergeClusters, "plotpicker-pzmap.gif", "Merge brushed cells to same group")
in_addbuttonpush(1, 0, @in_RenameCellTypeBatchID, "plotpicker-scatterhist.gif", "Rename cell type or batch ID");
in_addbuttonpush(1, 1, @in_ClusterCellsS, "plotpicker-dendrogram.gif", "Clustering using cell embedding (S)")
in_addbuttonpush(1, 0, @in_ClusterCellsX, "icon-mw-cluster-10.gif", "Clustering using expression matrix (X)")
in_addbuttonpush(1, 1, {@in_DetermineCellTypeClustersGeneral, true}, "plotpicker-contour.gif", "Assign cell types to groups")
in_addbuttonpush(1, 0, @in_Brush4Celltypes, "brush.gif", "Assign cell type to selected cells");
in_addbuttonpush(1, 1, @gui.callback_Brush4Markers, "plotpicker-kagi.gif", "Marker genes of brushed cells");
in_addbuttonpush(1, 0, @gui.callback_FindAllMarkers, "plotpicker-plotmatrix.gif", "Marker gene heatmap");
in_addbuttonpush(1, 1, @gui.callback_ShowClustersPop, "plotpicker-geoscatter.gif", "Show cell clusters/groups individually");
in_addbuttonpush(1, 0, @gui.callback_SelectCellsByClass, "plotpicker-pointfig.gif", "Select cells by class");
in_addbuttonpush(1, 0, @in_DeleteSelectedCells, "plotpicker-qqplot.gif", "Delete selected cells");
in_addbuttonpush(1, 0, @callback_SaveX, "export.gif", "Export & save data");
in_addbuttonpush(1, 1, @in_EmbeddingAgain, "plotpicker-geobubble.gif", "Embedding (tSNE, UMP, PHATE)");
in_addbuttonpush(1, 0, @in_Switch2D3D, "plotpicker-image.gif", "Switch 2D/3D");
in_addbuttonpush(1, 1, @callback_CloseAllOthers, "icon-fa-cut-10.gif", "Close all other figures");
in_addbuttonpush(1, 0, @callback_PickPlotMarker, "plotpicker-rose.gif", "Switch scatter plot marker type");
in_addbuttonpush(1, 0, @gui.callback_PickColorMap, "plotpicker-compass.gif", "Pick new color map");
in_addbuttonpush(1, 0, @in_RefreshAll, "icon-mat-refresh-20.gif", "Refresh");
in_addbuttonpush(0, 0, @gui.callback_MultiGroupingViewer, "plotpicker-arxtimeseries.gif", "Multi-grouping View...");
in_addbuttonpush(0, 0, @gui.callback_CrossTabulation, "plotpicker-comet.gif", "Cross tabulation");
in_addbuttonpush(0, 1, @gui.callback_Violinplot, "violinplot.gif", "Gene Violin Plot...");
in_addbuttonpush(0, 0, @gui.callback_DrawDotplot, "icon-mat-blur-linear-10.gif", "Gene Dot Plot...");
in_addbuttonpush(0, 0, @gui.callback_GeneHeatMap, "icon-mat-apps-20.gif", "Gene Heatmap...");
in_addbuttonpush(0, 1, @in_CompareGeneBtwCls, "cellscore2.gif", "Cell score analysis--obtaining gene signature score for each cell");
in_addbuttonpush(0, 0, @gui.callback_GetCellSignatureMatrix, "icon-fa-connectdevelop-20.gif", "Cell state analysis--obtaining multiple gene signature scores to reveal functional state of cells");
in_addbuttonpush(0, 1, @gui.callback_DEGene2Groups, "plotpicker-boxplot.gif", "Differential expression (DE) analysis)");
in_addbuttonpush(0, 0, @gui.callback_DPGene2Groups, "plotpicker_noisepsd.gif", "Differential program (DP) analysis)");
in_addbuttonpush(0, 0, @in_EnrichrHVGs, "plotpicker-andrewsplot.gif", "Functional enrichment analysis with HVGs");
in_addbuttonpush(0, 1, @gui.callback_BuildGeneNetwork, "noun_Network_691907.gif", "Build gene regulatory network");
in_addbuttonpush(0, 0, @gui.callback_CompareGeneNetwork, "noun_Deep_Learning_2424485.gif", "Compare two scGRNs");
in_addbuttonpush(0, 1, {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
gui.add_3dcamera(DeftToolbarHandle, 'AllCells');
drawnow;

m_fil = uimenu(FigureHandle,'Text','&File','Accelerator','F');
in_addmenu(m_fil, 0, @OpenSCEDataFilematMenuSelected, 'Open SCE Data File (*.mat)...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Open TXT/TSV/CSV File (*.txt)...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Open Seurat/Rds File (*.rds)...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Open AnnData/H5ad File (*.h5ad)...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Open Loom File (*.loom)...');

in_addmenu(m_fil, 1, @ExitMenuSelected, 'Open 10x Genomics H5 File (*.h5)...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Open 10x Genomics MTX File (*.mtx)...');
in_addmenu(m_fil, 1, @ExitMenuSelected, 'Import From 10x Genomics ''outs'' Folder...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Import From Parse Biosciences ''outs'' Folder...');
in_addmenu(m_fil, 1, @ExitMenuSelected, 'Read From Link to GEO h5 File...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Read From Link to GEO mtx.gz File...');
in_addmenu(m_fil, 0, @ExitMenuSelected, 'Read From Link to GEO txt.gz File...');
in_addmenu(m_fil, 0, @ReadFromGEOAccession, 'Read From GEO Accession Number(s)...');
in_addmenu(m_fil, 1, @LoadSCEFromWorkspace, 'Load SCE Variable from Workspace...');
in_addmenu(m_fil, 0, @LoadExampleDataMenuSelected, 'Load Example Data...');
in_addmenu(m_fil, 1, @CloseSCEDataMenuSelected, 'Close SCE Data');
in_addmenu(m_fil, 1, @ExitMenuSelected, 'Exit');

m_edi = uimenu(FigureHandle,'Text','&Edit','Accelerator','E');
in_addmenu(m_edi,0,@in_SelectCellsByQC,'Filter genes and cells');
in_addmenu(m_edi,1,@in_Brushed2NewCluster,'Add brushed cells to a new group');
in_addmenu(m_edi,0,@in_Brushed2MergeClusters,'Merge brushed cells to same group');
in_addmenu(m_edi,1,@in_RenameCellTypeBatchID,'Rename cell type or batch ID');
in_addmenu(m_edi,0,@in_RenameCellTypeBatchID,'Rename cell type or batch ID');
in_addmenu(m_edi, 1, @in_AddEditCellAttribs, 'Add/Edit Cell Attributes...');
in_addmenu(m_edi, 0, @in_ExportCellAttribTable, 'Export Cell Attribute Table...');
in_addmenu(m_edi, 0, @gui.callback_ViewMetaData, 'View Metadata...');
in_addmenu(m_edi, 1, {@in_MergeCellSubtypes, 1}, 'Import Cell Annotation from SCE in Workspace...');
in_addmenu(m_edi, 0, {@in_MergeCellSubtypes, 2}, 'Import Cell Annotation from SCE Data File...');

m_vie = uimenu(FigureHandle,'Text','&View','Accelerator','V');
in_addmenu(m_vie,0,@gui.callback_ShowGeneExpr,'Select genes to show expression');
in_addmenu(m_vie,0,@in_ShowCellStates,'Show cell state');
in_addmenu(m_vie,0,@in_SelectCellsByQC,'Label cell groups');
in_addmenu(m_vie,1,@gui.callback_MultiEmbeddingViewer,'Multi-embedding View...');
in_addmenu(m_vie,0,@gui.callback_MultiGroupingViewer,'Multi-grouping View...');
in_addmenu(m_vie,0,@gui.callback_CrossTabulation,'Cross Tabulation...');
in_addmenu(m_vie, 0, @gui.callback_ShowHgBGeneExpression, 'Show Hemoglobin (Hgb) Genes Expression...');
in_addmenu(m_vie, 0, @gui.callback_ShowMtGeneExpression, 'Show Mitochondrial (Mt-) Genes Expression...');

m_net = uimenu(FigureHandle, 'Text', '&Network', 'Accelerator', 'N');
in_addmenu(m_net, 0, @in_Select5000Genes, 'Remove Less Informative Genes to Reduce Gene Space...');
in_addmenu(m_net, 0, @gui.i_setnetwd, 'Set Network Analysis Working Root Directory...');
in_addmenu(m_net, 1, {@in_scTenifoldNet,1}, 'Construct GRN using PC Regression [PMID:33336197] 🐢...');
in_addmenu(m_net, 0, {@in_scTenifoldNet,2}, 'Construct & Compare GRNs (scTenifoldNet Analysis) [PMID:33336197] 🐢...');
in_addmenu(m_net, 1, @gui.callback_scTenifoldKnk1, 'Virtual Gene KO - scTenifoldKnk [PMID:35510185] 🐢🐢 ...');
in_addmenu(m_net, 0, @gui.callback_VirtualKOGenKI, 'Virtual Gene KO - GenKI [PMID:37246643] (Python Required) 🐢🐢 ...');
in_addmenu(m_net, 1, @gui.callback_scTenifoldXct, 'Cell-Cell Interactions (CCIs) - scTenifoldXct [PMID:36787742] 🐢🐢 ...');
in_addmenu(m_net, 0, @gui.callback_scTenifoldXct2, 'Differential CCIs - scTenifoldXct [PMID:36787742] 🐢🐢🐢 ...');

m_ext = uimenu(FigureHandle, 'Text', 'E&xternal', 'Accelerator', 'x');
in_addmenu(m_ext, 0, @gui.i_setrenv, 'Check R Environment');
in_addmenu(m_ext, 0, @gui.i_setpyenv, 'Check Python Environment');
in_addmenu(m_ext, 0, @gui.i_setextwd, 'Set External Program Working Root Directory...');
in_addmenu(m_ext, 1, @in_DecontX, 'Detect Ambient RNA Contamination (DecontX/R) [PMID:32138770]...');
in_addmenu(m_ext, 0, @in_RunSeuratWorkflow, 'Run Seurat/R Workflow (Seurat/R) [PMID:25867923]...');
in_addmenu(m_ext, 0, @gui.callback_RunMonocle3, 'Pseudotime Analysis (Monocle3/R) [PMID:28825705]...');
in_addmenu(m_ext, 1, @gui.callback_MELDPerturbationScore, 'MELD Perturbation Score (MELD/Py) [PMID:33558698]...');
in_addmenu(m_ext, 0, {@in_SubsampleCells, 2}, 'Geometric Sketching (geosketch/Py) [PMID:31176620]...');
in_addmenu(m_ext, 0, @in_HarmonyPy, 'Batch Integration (Harmony/Py) [PMID:31740819]...');
in_addmenu(m_ext, 0, @in_DoubletDetection, 'Detect Doublets (Scrublet/Py) [PMID:30954476]...');
in_addmenu(m_ext, 0, @in_RunDataMapPlot, 'Run DataMapPlot (datamapplot/Py)...');
in_addmenu(m_ext, 1, @gui.callback_ExploreCellularCrosstalk, 'Talklr Intercellular Crosstalk [DOI:10.1101/2020.02.01.930602]...');

m_tol = uimenu(FigureHandle, 'Text', '&Tools', 'Accelerator', 'T');
in_addmenu(m_tol, 0, @gui.callback_SelectCellsByMarker, 'Extract Cells by Marker (+/-) Expression...');
in_addmenu(m_tol, 0, @in_MergeSubCellTypes, 'Merge Subclusters of the Same Cell Type');
in_addmenu(m_tol, 1, @in_WorkonSelectedGenes, 'Select Top n  Highly Variable Genes (HVGs) to Work on...');
in_addmenu(m_tol, 0, @in_SubsampleCells, 'Subsample 50% Cells to Work on...');
in_addmenu(m_tol, 1, @gui.callback_DEGene2GroupsBatch, 'Differential Expression (DE) Analysis in Cell Type Batch Mode...');
in_addmenu(m_tol, 0, @gui.callback_DPGene2GroupsBatch, 'Differential Program (DP) Analysis in Cell Type Batch Mode...');
in_addmenu(m_tol, 1, @gui.callback_CalculateGeneStats, 'Calculate Gene Expression Statistics...');
in_addmenu(m_tol, 0, @gui.callback_CellCycleLibrarySize, 'Library Size of Cell Cycle Phases...');
in_addmenu(m_tol, 0, @gui.callback_CellCycleAssignment, 'Assign Cell Cycle Phase...');
in_addmenu(m_tol, 1, {@in_DetermineCellTypeClustersGeneral, false}, 'Annotate Cell Type Using Customized Markers...');
in_addmenu(m_tol, 0, @in_SubtypeAnnotation, 'Annotate Cell Subtype...');

m_exp = uimenu(FigureHandle, 'Text', 'Ex&perimental', 'Accelerator', 'p');
in_addmenu(m_exp, 1, @gui.callback_SplitAtacGex, 'Split Multiome ATAC+GEX Matrix...');
in_addmenu(m_exp, 0, @in_DrawKNNNetwork, 'Plot Cell kNN Network...');
in_addmenu(m_exp, 0, @in_DrawTrajectory, 'Plot Cell Trajectory...');
in_addmenu(m_exp,0,{@MergeCellSubtypes,1,true},'Import All Cell Annotation from SCE in Workspace...');
in_addmenu(m_exp,0,{@MergeCellSubtypes,2,true},'Import All Cell Annotation from SCE Data File...');
in_addmenu(m_exp, 1, {@in_MergeSCEs, 1}, 'Merge SCE Variables in Workspace...');
in_addmenu(m_exp, 0, {@in_MergeSCEs, 2}, 'Merge SCE Data Files...');
in_addmenu(m_exp, 0, {@gui.i_savemainfig, 2}, 'Save Figure as Graphic File...');
in_addmenu(m_exp, 0, {@gui.i_savemainfig, 1}, 'Save Figure as SVG File...');
in_addmenu(m_exp, 1, @in_SingleClickSolution, 'Single Click Solution (from Raw Data to Annotation)...');

m_hlp = uimenu(FigureHandle, 'Text', '&Help', 'Accelerator', 'H');
in_addmenu(m_hlp, 1, {@(~, ~) web('https://scgeatool.github.io/')}, 'Visit SCGEATOOL-Standalone Website...');
in_addmenu(m_hlp, 0, @callback_CheckUpdates, 'Check for Updates...');
drawnow;

% handles = guihandles( FigureHandle );
guidata( FigureHandle, sce );
FigureHandle.Visible="on";
in_update_figure;
drawnow;

if useexample
    close(pD);
end





    function [out]=in_gscatter3(fig, s,c)
        if size(s,2)>=3
            [out]=scatter3(fig, s(:,1), s(:,2), s(:,3), 10, c);
        else
            [out]=scatter(fig, s(:,1), s(:,2), 10, c);
        end
    end

    function in_update_figure
        focus(FigureHandle);
        if sce.NumCells>0
            % kc = numel(unique(c));
            % colormap(pkg.i_mycolorlines(kc));
            % c=sce.c;
            %[ax,bx]=view(hAx);
            mean(sce.c)
            [h]=in_gscatter3(hAx,sce.s,sce.c);
            %view(hAx,ax,bx);
            title(hAx, sce.title);
            subtitle(hAx, '[genes x cells]');
            grid(hAx,"on");
            drawnow;
            hAx.Visible="on";
        else
            if ~isempty(h) && isvalid(h)
                h.Visible="off";
            end
            hAx.Visible="off";
        end
    end

    function ExitMenuSelected(~, ~)
        delete(FigureHandle)
    end

    function CloseSCEDataMenuSelected(~, ~)
        sce = SingleCellExperiment;
        %delete(hAx.Children);
        %hAx.Visible="off";
        in_update_figure;
    end

    function OpenSCEDataFilematMenuSelected(~, ~)
        [fname, pathname] = uigetfile( ...
                        {'*.mat', 'SCE Data Files (*.mat)'; ...
                        '*.*', 'All Files (*.*)'}, ...
                        'Pick a SCE Data File');
        if isequal(fname, 0), return; end
        scefile = fullfile(pathname, fname);
        try                
            load(scefile, 'sce');            
            in_update_figure;
        catch ME            
            uialert(FigureHandle,ME.message,"");
            return;
        end
    end

    function ReadFromGEOAccession(~,~)
        acc = inputdlg({'Input Number(s) (e.g., GSM3308547,GSM3308548):'}, ...
            'GEO Accession', [1, 50], {'GSM3308547'});
        if isempty(acc), return; end
        %acc = strtrim(deblank(acc{1}));
        %acc = strrep(acc,' ','');
        acc = regexprep(acc{1},'[^a-zA-Z0-9,;]','');
        if isempty(acc) || ~strlength(acc) > 4, return; end
        if strlength(acc) > 4 && ~isempty(regexp(acc, 'G.+', 'once'))
            accv = unique(strsplit(acc, {',', ';', ' '}), 'stable');
            if length(accv) > 1
                dmanswer = questdlg('Download and merge data sets?', ...
                    '', 'Yes', 'Cancel', 'Yes');
                if ~strcmp(dmanswer, 'Yes'), return; end
                try
                    fw = uiwaitbar(FigureHandle);
                    [sce] = pkg.pipeline_multisamplesmerge(accv, false);
                    uiwaitbar(FigureHandle,fw);
                catch ME
                    uiwaitbar(FigureHandle,fw,true);
                    uialert(FigureHandle, ME.message,"");
                    return;
                end
            else                        
                try
                    fw = uiwaitbar(FigureHandle);
                    [sce] = sc_readgeoaccession(acc);
                    uiwaitbar(FigureHandle,fw);
                catch ME
                    uiwaitbar(FigureHandle,fw,true);
                    uialert(FigureHandle, ME.message,"");
                    return;
                end
            end
            %                     metainfo=sprintf("Source: %s",acc);
            %                     sce=sce.appendmetainfo(metainfo);
        end
        if isempty(sce), return; end
        if length(sce.g) ~= length(unique(sce.g))
            disp('Construct unique gene names from input gene list.')
            sce.g = matlab.lang.makeUniqueStrings(sce.g);
        end
        in_update_figure;
    end

    function LoadSCEFromWorkspace(~,~)
        a = evalin('base', 'whos');
        b = struct2cell(a);
        valididx = ismember(b(4, :), 'SingleCellExperiment');
        if isempty(valididx) || ~any(valididx)
            uialert(FigureHandle,'No SCE in the Workspace.', '');
            return;
        end
        a = a(valididx);
        [indx, tf] = listdlg('PromptString', {'Select SCE variable:'}, ...
            'liststring', b(1, valididx), 'SelectionMode', 'multiple');
        if tf == 1
            if length(indx) == 1
                sce = evalin('base', a(indx).name);
            elseif length(indx) > 1
                answer = questdlg('Which set operation method to merge genes?', 'Merging method', ...
                    'Intersect', 'Union', 'Intersect');
                if ~ismember(answer, {'Union', 'Intersect'}), return; end
                methodtag = lower(answer);
                try
                    insce = cell(1, length(indx));
                    s = "";
                    for k = 1:length(indx)
                        insce{k} = evalin('base', a(indx(k)).name);
                        s = sprintf('%s,%s', s, a(indx(k)).name);
                    end
                    s = s(2:end);
                    fprintf('>> sce=sc_mergesces({%s},''%s'');\n', s, methodtag);
                    fw = uiwaitbar(FigureHandle);
                    sce = sc_mergesces(insce, methodtag);
                catch ME
                    uiwaitbar(FigureHandle,fw);
                    uialert(FigureHandle, ME.message, "");
                    return;
                end
                gui.gui_waitbar(fw);
            end
        else
            return;
        end
        in_update_figure;
    end

    function LoadExampleDataMenuSelected(~, ~)        
        selection = uiconfirm(FigureHandle, ...
            "Load processed or raw data?", "Load Data", ...
            "Options",["Processed","Raw","Cancel"], ...
            "DefaultOption",1,"CancelOption",3);
        switch selection 
            case 'Processed'
                answer=1;
            case 'Raw'
                answer=2;
            case 'Cancel'
                return;
        end
        pw1 = fileparts(mfilename('fullpath'));
        fprintf('Loading SCE Data File example_data/workshop_example.mat...');
        tic;
        file1 = fullfile(pw1, '..', 'example_data', 'workshop_example.mat');
        if ~exist(file1, "file")
            uialert(FigureHandle,"Example data file does not exist.","Missing File");
            return;
        end
        load(file1, 'sce');
        if answer == 2
            a = SingleCellExperiment(sce.X, sce.g);
            a.c_batch_id = sce.c_batch_id;
            a.c_cell_id = sce.c_cell_id;
            a.metadata = sce.metadata;
            sce = a;
        end
        fprintf('Done.\n');
        toc;
        in_update_figure;
    end



% xxxxxxx

    function in_callback_ShowGeneExpr(~,~)
        [selectedg] = ui_dropdown(FigureHandle,sce.g);
        if ~empty(selectedg)
            sce.c = sce.X(sce.g == selectedg, :);
            in_update_figure;
        end
    end


    function in_fixfield(oldf,newf)
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

        if askpref
            %  gui.gui_userguidingpref(false);
            %answer=uiconfirm(FigureHandle, 'Show User Onboarding Toolbar again next time?','');
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
        end
        pt = uipushtool(barhandle, 'Separator', sepTag);
        %pt.Icon = in_getPtImage(imgFil);
        pt.Icon = fullfile(mfolder,'..','resources',imgFil);
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
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
        end
        pt = uitoggletool(barhandle, 'Separator', sepTag);        
        pt.Icon = fullfile(mfolder,'..','resources',imgFil);
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
    end

    function in_togglebtfun(src, ~, func, ~, imgFil, ...
            actiondelay, tooltipTxt)
        if nargin < 6, actiondelay = true; end
        src.Icon = fullfile(mfolder,'..','resources',imgFil);
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

    function in_closeRequest(hObject, ~)
        if ~(ismcc || isdeployed)
            ButtonName = uiconfirm(FigureHandle, 'Save SCE before closing SCGEATOOL?','');
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

    function in_GEOAccessionToSCE(src, ~)
        answer = uiconfirm(FigureHandle, 'Current SCE will be replaced. Continue?','');
        if ~strcmp(answer, 'Yes'), return; end       
        acc = inputdlg({'Input Number(s) (e.g., GSM3308547,GSM3308548):'}, ...
                    'GEO Accession', [1, 50], {'GSM3308547'});        
        if isempty(acc), return; end
        acc = acc{1};
        if strlength(acc) > 4 && ~isempty(regexp(acc, 'G.+', 'once'))
            try
                fw = uiwaitbar(FigureHandle);
                [sce] = sc_readgeoaccession(acc);
                [c, cL] = grp2idx(sce.c);
                uiwaitbar(FigureHandle, fw);
                guidata(FigureHandle, sce);
                in_RefreshAll(src, [], false, false);
            catch ME
                uiwaitbar(FigureHandle, fw);
                uialert(FigureHandle, ME.message, "");
            end
        end
    end

    function in_RunDataMapPlot(src, ~)        
        ndim = 2;
        [vslist] = gui.i_checkexistingembed(sce, ndim);
        if isempty(h.ZData) && size(sce.s,2)==2 && length(vslist) <= 1
            gui.callback_RunDataMapPlot(src, []);
        elseif isempty(h.ZData) && size(sce.s,2)==2 && length(vslist) > 1
            answer = uiconfirm(FigureHandle, 'Using current 2D embedding?','');
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
            answer=uiconfirm(FigureHandle, ...
                'This function requires 2D embedding. Continue?','');
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
            answer = uiconfirm(FigureHandle, 'Import annotation for all cells or just cells of a subtype?', '', ...
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
            k = gui.i_inputnumk(2000, 1, sce.NumGenes, 'the number of HVGs');
            if isempty(k), return; end
            answer = uiconfirm(FigureHandle, 'Which method?', 'Select Method', ...
                'Brennecke et al. (2013)', 'Splinefit Method', ...
                'Brennecke et al. (2013)');
            fw = uiwaitbar(FigureHandle);
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
            if ~all(y), uialert(FigureHandle, 'Runtime error.',"");
                return;
            end
            sce.g = sce.g(idx);
            sce.X = sce.X(idx, :);
            uiwaitbar(FigureHandle, fw);
            in_RefreshAll(src, [], true, false);
    end

    function in_SubsampleCells(src, ~, methodoption)
        if nargin < 3
            methodoption = [];
        end
        answer1 = uiconfirm(FigureHandle, 'This function subsamples 50% of cells. Continue?','');
        if ~strcmp(answer1, 'OK')
            return;
        end

        if isempty(methodoption)
            answer = uiconfirm(FigureHandle, 'Select method:', '', ...
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
            fw = uiwaitbar(FigureHandle);
            Xn = log(1+sc_norm(sce.X))';
            [~, Xn] = pca(Xn, 'NumComponents', 300);
            uiwaitbar(FigureHandle, fw);
            try
                ids = run.py_geosketch(Xn, tn);
            catch ME
                uiwaitbar(FigureHandle, fw, true);
                uialert(FigureHandle, ME.message,"");
                return;
            end

        end
        if ~isempty(ids)
            sce = sce.selectcells(ids);
            c = sce.c;
            in_RefreshAll(src, [], true, false);
        else
            uialert(FigureHandle, 'Running error. No action is taken.',"");
        end
    end

    function in_SingleClickSolution(src, ~)
        speciestag = gui.i_selectspecies(2);
        if isempty(speciestag), return; end

        fw = uiwaitbar_adv;
        uiwaitbar_adv(fw,1/8,'Basic QC Filtering...');
        sce = sce.qcfilter;
        uiwaitbar_adv(fw,2/8, 'Embeding Cells Using tSNE...');
        sce = sce.embedcells('tsne3d', true, true, 3);

        uiwaitbar_adv(fw,3/8, 'Clustering Cells Using K-means...');
        sce = sce.clustercells([], [], true);
        uiwaitbar_adv(fw,4/8, 'Annotating Cell Type Using PanglaoDB...');
        
        tic
        sce = sce.assigncelltype(speciestag, false);
        toc
        uiwaitbar_adv(fw,5/8, 'Estimate Cell Cycles...');
        
        sce = sce.estimatecellcycle;
        uiwaitbar_adv(fw,6/8, 'Estimate Differentiation Potency of Cells...');

        sce = sce.estimatepotency(speciestag);

        uiwaitbar_adv(fw,7/8);
        [c,cL] = grp2idx(sce.c_cell_type_tx);
        sce.c = c;
        uiwaitbar_adv(fw);
        in_RefreshAll(src, [], true, false);
        ix_labelclusters(true);
        setappdata(FigureHandle, 'cL', cL);
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
            uialert(FigureHandle, ME.message,"");
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
        fw = uiwaitbar(FigureHandle);
        try
            [sce] = run.r_seurat(sce, ndim, wkdir);
            [c, cL] = grp2idx(sce.c);
        catch
            uiwaitbar(FigureHandle, fw);
            return;
        end
        uiwaitbar(FigureHandle, fw);
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
        fw = uiwaitbar(FigureHandle);
        try
            [Xdecon, contamination] = run.r_decontX(sce, wkdir);
        catch
            uiwaitbar(FigureHandle, fw, true);
            uialert(FigureHandle, 'Runtime error.',"");
            return;
        end
        uiwaitbar(FigureHandle, fw);
        figure('WindowStyle', 'modal');
        gui.i_stemscatter(sce.s, contamination);
        % zlim([0 1]);
        zlabel('Contamination rate')
        title('Ambient RNA contamination')
        answer = uiconfirm(FigureHandle, "Remove contamination?",'');
        switch answer
            case 'OK'
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
            answer = uiconfirm(FigureHandle, 'Color cells by batch id (SCE.C_BATCH_ID)?', '');
            switch answer
                case 'OK'
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
    
            ButtonName = uiconfirm(FigureHandle, 'Update Saved Embedding?', '');
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
            answer = uiconfirm(FigureHandle, sprintf("Remove %d doublets?", sum(isDoublet)),'');
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

    function in_RefreshAll(src, ~, keepview, keepcolr)
        if nargin < 4, keepcolr = false; end
        if nargin < 3, keepview = false; end
        if keepview || keepcolr
            [para] = gui.i_getoldsettings(src);
        end
        figure(FigureHandle);
        % was3d = ~isempty(h.ZData);
        if size(sce.s, 2) >= 3
            if keepview, [ax, bx] = view(hAx); end
            [h] = in_gscatter3(hAx, sce.s, c);
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
            kc = numel(unique(c));
            colormap(pkg.i_mycolorlines(kc));
        end
        title(hAx, sce.title);
        subtitle(hAx, '[genes x cells]');
    end

    function in_Switch2D3D(src, ~)  
        [para] = gui.i_getoldsettings(src);

        if isempty(h.ZData)               % current 2D xxx
            ansx = uiconfirm(FigureHandle, 'Switch to 3D?','');
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
                    ansx = uiconfirm(FigureHandle, ...
                        'Using existing 3D embedding? Select "No" to re-embed.','');
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
            ansx = uiconfirm(FigureHandle, 'Switch to 2D?','');
            if ~strcmp(ansx, 'Yes'), return; end
            [vslist] = gui.i_checkexistingembed(sce, 2);
            if ~isempty(vslist)
                answer = uiconfirm(FigureHandle, 'How to make 2D embedding?','', ...
                    'Pick existing 2D','Re-embed cells', ...
                    'Reduce current 3D','Pick existing 2D');
            else
                answer = uiconfirm(FigureHandle, 'How to make 2D embedding?','', ...
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
                    answer2 = uiconfirm(FigureHandle, 'Which view to be used to project cells?', '', ...
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
        answer = uiconfirm(FigureHandle, 'Add or edit cell attribute?','', ...
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
            answer = uiconfirm(FigureHandle, ...
                'Rename cell type, batch ID, or gene name?', ...
                '', "Options",["Cell type","Batch ID","Gene name"], ...
                "DefaultOption",2,"CancelOption",3);
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
            answer = uiconfirm(FigureHandle, 'Using exsiting embedding?','');
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
            %     answer1 = uiconfirm(FigureHandle, sprintf('Use existing %s embedding or re-compute new embedding?', ...
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
                fw = uiwaitbar(FigureHandle);
                try
                    forced = true;
                    %if contains(methoddimtag, 'tsne'), disp('tSNE perplexity = 30'); end
                    sce = sce.embedcells(methodtag, forced, usehvgs, ndim, K);
                    % disp('Following the library-size normalization and log1p-transformation, we visualized similarity among cells by projecting them into a reduced dimensional space using t-distributed stochastic neighbor embedding (t-SNE)/uniform manifold approximation and projection (UMAP).')
                catch ME
                    uiwaitbar(FigureHandle,fw, true);
                    uialert(FigureHandle, ME.message,"");
                    % rethrow(ME)
                    return;
                end
                uiwaitbar(FigureHandle, fw);
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
        if ~manuallyselect, fw = uiwaitbar_adv; end

        for ix = 1:max(c)
            if ~manuallyselect
                uiwaitbar_adv(fw, ix/max(c));
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
        if ~manuallyselect, uiwaitbar_adv(fw); end
        sce.c_cell_type_tx = string(cL(c));
        nx = length(unique(sce.c_cell_type_tx));
        if nx > 1
            newtx = erase(sce.c_cell_type_tx, "_{"+digitsPattern+"}");
            if length(unique(newtx)) ~= nx
                answer = uiconfirm(FigureHandle, ...
                    'Merge subclusters of same cell type?','');
                if strcmp(answer, 'Yes')
                    in_MergeSubCellTypes(src);
                end
            end
        end
        guidata(FigureHandle, sce);
    end

    function [iscelltype] = in_pickcelltypeclusterid(a)
        answer = uiconfirm(FigureHandle, a, '', ...
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
        answer = uiconfirm(FigureHandle, ...
            'This is a one-time analysis. Cell type labels will not be saved. Continue?','');
        if ~strcmp(answer, 'Yes')
            return;
        end
        speciestag = gui.i_selectspecies;
        if isempty(speciestag), return; end
        fw = uiwaitbar(FigureHandle);
        [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, sce.s, ...
            ptsSelected, ...
            speciestag, "all", "panglaodb", false);
        ctxt = Tct.C1_Cell_Type;
        uiwaitbar(FigureHandle, fw);

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
        answer1 = uiconfirm(FigureHandle, 'Display in place or in new figure?', '', ...
            'Options',["In place", "New figure","Cancel"]);            
            
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
        
        answer = uiconfirm(FigureHandle, ...
            'This function applies to a homogeneous group of cells. Remove lowly expressed genes before applying. Continue?','');
        if ~strcmp(answer, 'OK'), return; end
        
        ptsSelected = logical(h.BrushData.');
        if any(ptsSelected)
            [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, sce);
            if ~letdoit, return; end
            if sum(ptsSelected) < 200
                answer = uiconfirm(FigureHandle, sprintf('Too few cells (n = %d) selected, continue?', sum(ptsSelected)));
                if ~strcmp(answer, 'Yes'), return; end
            end           
            scetmp = sce.removecells(~ptsSelected);
            scetmp = scetmp.qcfilter(1000, 0.15, 15);
            gui.callback_EnrichrHVGs(src, events, scetmp);
        else
            answer = uiconfirm(FigureHandle, sprintf('All cells (n = %d) included, continue?', sce.NumCells));
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

        answer = uiconfirm(FigureHandle, 'Delete selected or unselected cells?', '', ...
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
            fw = uiwaitbar(FigureHandle);
        end
        sce = sce.removecells(ptsSelected);
        if needprogressbar
            uiwaitbar(FigureHandle, fw);
        end
        [c, cL] = grp2idx(sce.c);
        in_RefreshAll(src, [], true, true);
    end

    function in_scTenifoldNet(src,events,methodtag)
        if numel(unique(sce.c_cell_type_tx))>1
            answer=uiconfirm(FigureHandle, 'This analysis is cell type-specific; however, current SCE contains multiple cell types. Continue?');
            if ~strcmp(answer,'Yes'), return; end
        end
        answer=uiconfirm(FigureHandle, 'Subsample cells?','','Yes 🐢','No 🐇','Cancel','No 🐇');
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
        fw = uiwaitbar(FigureHandle);
        set(0, 'CurrentFigure', FigureHandle);
        figure('WindowStyle', 'modal');
        sc_knngraph(sce.s, k, true);
        uiwaitbar(FigureHandle, fw);
    end

    function in_DrawTrajectory(src, ~)
        % waitfor(warndlg('This function should not be applied to tSNE and UMAP embeddings, as they "encourage a representation of the data as disjoint clusters, which is less meaningful for modeling continuous developmental trajectories" [PMID:25664528].', ''));
        if ~isempty(f_traj)
            answer = uiconfirm(FigureHandle, 'Remove existing trajectory curve?');
            switch answer
                case 'Yes'
                    in_RefreshAll(src, [], true, true);  % keepview, keepcolr
                case 'No'
                otherwise
                    return;
            end
        end
        answer = uiconfirm(FigureHandle, 'Which method?', '', 'splinefit', 'princurve', ...
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
        % answer = uiconfirm(FigureHandle, 'Save/Update pseudotime T in SCE', ...
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

        answer = uiconfirm(FigureHandle, 'View expression of selected genes', ...
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
        answer = uiconfirm(FigureHandle, 'Cluster cells using embedding S?');
        if ~strcmp(answer, 'Yes'), return; end

        [sx] = gui.i_pickembedvalues(sce);
        if isempty(sx), return; end

        answer = uiconfirm(FigureHandle, 'Which method?', 'Select Algorithm', ...
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
        answer = uiconfirm(FigureHandle, 'Cluster cells using expression matrix X?');
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
            answer1 = uiconfirm(FigureHandle, sprintf('Using existing %s clustering?', ...
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
            fw = uiwaitbar(FigureHandle);
            try
                % [sce.c_cluster_id]=sc_cluster_x(sce.X,k,'type',methodtag);
                sce = sce.clustercells(k, methodtag, true, sx);
            catch ME
                uiwaitbar(FigureHandle, fw, true);
                uialert(FigureHandle, ME.message,"");
                return
            end
            uiwaitbar(FigureHandle, fw);
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
            sce = guidata(FigureHandle);
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
            setappdata(FigureHandle, 'cL', cL);
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
                    answer = uiconfirm(FigureHandle, sprintf('Label %d groups with index or text?', ...
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