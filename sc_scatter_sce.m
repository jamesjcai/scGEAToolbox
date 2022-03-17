function varargout = sc_scatter_sce(sce, varargin)

if usejava('jvm') && ~feature('ShowFigureWindows')
    error('MATLAB is in a text mode. This function requires a GUI-mode.');
end
if nargin < 1
    % error('Usage: sc_scatter_sce(sce)');
    sc_scatter;
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
parse(p, sce, varargin{:});
cin = p.Results.c;
sin = p.Results.s;
methodid = p.Results.methodid;
ax = [];
bx = [];
tmpcelltypev=cell(sce.NumCells,1);

if isempty(cin)
    sce.c = ones(size(sce.X, 2), 1);
else
    sce.c = cin;
end
if ~isempty(sin)
    sce.s = sin;
end

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

hAx = axes('Parent', FigureHandle);

[h] = gui.i_gscatter3(sce.s, c, methodid,1,hAx);
title(hAx,sce.title);

dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};



defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle

UitoolbarHandle = uitoolbar('Parent', FigureHandle);
set(UitoolbarHandle, 'Tag', 'FigureToolBar', ...
    'HandleVisibility', 'off', 'Visible', 'on');

mfolder = fileparts(mfilename('fullpath'));


i_addbutton(1,0,@callback_ShowGeneExpr,"list.gif","Select a gene to show expression")
i_addbutton(1,0,@ShowCellStates,"list2.gif","Show cell state")
i_addbutton(1,0,@SelectCellsByQC,"plotpicker-effects.gif","Filter genes and cells")

%i_addbutton(1,1,@LabelClusters,"plotpicker-scatter.gif","Label clusters")
ptlabelclusters = uitoggletool(UitoolbarHandle, 'Separator', 'on');
[img, map] = imread(fullfile(mfolder, 'resources', 'plotpicker-scatter.gif'));
%[img, map] = imread(fullfile(matlabroot, ...
%    'toolbox', 'matlab', 'icons', 'plotpicker-scatter.gif'));
% map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
ptImage = ind2rgb(img, map);
ptlabelclusters.CData = ptImage;
ptlabelclusters.Tooltip = 'Label clusters';
ptlabelclusters.ClickedCallback = @LabelClusters;

i_addbutton(1,0,@Brushed2NewCluster,"plotpicker-glyplot-face.gif","Add brushed cells to a new cluster")
i_addbutton(1,0,@Brushed2MergeClusters,"plotpicker-pzmap.gif","Merge brushed cells to same cluster")
i_addbutton(1,1,@ClusterCellsS,"plotpicker-dendrogram.gif","Clustering using embedding S")
i_addbutton(1,0,@ClusterCellsX,"plotpicker-gscatter.gif","Clustering using expression matrix X")
i_addbutton(1,1,@DetermineCellTypeClusters,"plotpicker-contour.gif","Cell types of clusters")
i_addbutton(1,0,@Brush4Celltypes,"brush.gif","Cell types of brushed cells")
i_addbutton(1,0,@RenameCellType,"plotpicker-scatterhist.gif","Rename cell type")
i_addbutton(1,0,@callback_CellTypeMarkerScores,"cellscore.gif","Calculate Cell Scores from Cell Type Markers")
i_addbutton(1,0,@ShowCellStemScatter,"IMG00067.GIF","Stem scatter plot")
i_addbutton(1,1,@callback_Brush4Markers,"plotpicker-kagi.gif","Marker genes of brushed cells")
i_addbutton(1,0,@callback_MarkerGeneHeatmap,"plotpicker-plotmatrix.gif","Marker gene heatmap")
i_addbutton(1,1,@callback_ShowClustersPop,"plotpicker-geoscatter.gif","Show cell clusters/groups individually")
i_addbutton(1,0,@callback_SelectCellsByClass,"plotpicker-pointfig.gif","Select cells by class")
i_addbutton(1,0,@DeleteSelectedCells,"plotpicker-qqplot.gif","Delete selected cells")
i_addbutton(1,0,@callback_SaveX,"export.gif","Export & save data")
i_addbutton(1,1,@EmbeddingAgain,"plotpicker-geobubble.gif","Embedding (tSNE, UMP, PHATE)")
i_addbutton(1,0,@Switch2D3D,"plotpicker-image.gif","Switch 2D/3D")
i_addbutton(1,1,@callback_CloseAllOthers,"noun_Pruners_2469297.gif","Close All Other Figures")
i_addbutton(1,0,@callback_PickPlotMarker,"plotpicker-rose.gif","Switch scatter plot marker type")
i_addbutton(1,0,@callback_PickColorMap,"plotpicker-compass.gif","Switch color maps")
i_addbutton(1,0,@RefreshAll,"plotpicker-geobubble2.gif","Refresh")

i_addbutton(2,0,@call_scgeatool,"IMG00107.GIF"," ")
i_addbutton(2,0,@callback_CalculateCellScores,"cellscore2.gif","Calculate Cell Scores from List of Feature Genes")
i_addbutton(2,0,@callback_ComparePotency,"plotpicker-candle.gif","Compare Differentiation Potency");

i_addbutton(2,1,@gui.callback_MultiGroupingViewer,"plotpicker-arxtimeseries.gif","Multi-grouping View...");
i_addbutton(2,0,@gui.callback_CrossTabulation,"plotpicker-comet.gif","Cross Tabulation");

i_addbutton(2,1,@callback_CompareGeneBtwCls,"plotpicker-priceandvol.gif","Compare Gene Expression between Classes");
i_addbutton(2,0,@callback_DEGene2Groups,"plotpicker-boxplot.gif","Compare 2 groups (DE analysis)");
i_addbutton(2,0,@callback_GSEA_HVGs,"plotpicker-andrewsplot.gif","Function enrichment of HVG genes");
i_addbutton(2,0,@callback_BuildGeneNetwork,"noun_Network_691907.gif","Build gene regulatory network");
i_addbutton(2,0,@callback_CompareGeneNetwork,"noun_Deep_Learning_2424485.gif","Compare two scGRNs");

gui.add_3dcamera(defaultToolbar, 'AllCells');

m_vie = uimenu(FigureHandle,'Text','Multiv&iew');
m_vie.Accelerator = 'i';
i_addmenu(m_vie,0,@gui.callback_MultiEmbeddingViewer,'Multi-embedding View...');
i_addmenu(m_vie,0,@gui.callback_MultiGroupingViewer,'Multi-grouping View...');
i_addmenu(m_vie,0,@gui.callback_CrossTabulation,'Cross Tabulation...');
%uimenu(m_vie,'Text','Multi-embedding View...',...
%    'Callback',@gui.callback_MultiEmbeddingViewer);
%uimenu(m_vie,'Text','Multi-grouping View...',...    
%    'Callback',@gui.callback_MultiGroupingViewer);
%uimenu(m_vie,'Text','Cross Tabulation...',...
%    'Callback',@callback_CrossTabulation);

m_ext = uimenu(FigureHandle,'Text','Exte&rnal');
m_ext.Accelerator = 'r';
i_addmenu(m_ext,0,@gui.i_setrenv,'Check R Environment');
i_addmenu(m_ext,0,@gui.i_setpyenv,'Check Python Environment');
i_addmenu(m_ext,1,@DecontX,'Detect Ambient RNA Contamination (decontX/R required)...');

%uimenu(m_ext,'Text','Check R Environment',...
%    'Callback',@gui.i_setrenv);
%uimenu(m_ext,'Text','Check Python Environment',...
%    'Callback',@gui.i_setpyenv);
%uimenu(m_ext,'Text','Detect Ambient RNA Contamination (decontX/R required)...',...
%    'Separator','on',...        
%    'Callback',@DecontX);
i_addmenu(m_ext,0,@callback_SingleRCellType,'SingleR Cell Type Annotation (SingleR/R required)...');
i_addmenu(m_ext,0,@callback_RevelioCellCycle,'Revelio Cell Cycle Analysis (Revelio/R required)...');
i_addmenu(m_ext,0,@RunSeuratWorkflow,'Run Seurat/R Workflow (Seurat/R required)...');
i_addmenu(m_ext,0,@callback_MELDPerturbationScore,'MELD Perturbation Score (MELD/Python required)...');
i_addmenu(m_ext,0,@HarmonyPy,'Batch Integration (Harmony/Python required)...');
i_addmenu(m_ext,0,@DoubletDetection,'Detect Doublets (Scrublet/Python required)...');

m_exp = uimenu(FigureHandle,'Text','E&xperimental');
m_exp.Accelerator = 'x';
m_exp2 = uimenu(m_exp,'Text','sc&Tenifold Suite','Accelerator','T');

i_addmenu(m_exp2,0,@callback_scTenifoldNet1,'scTenifoldNet Construction 🐢🐢 ...');

%uimenu(m_exp2,'Text','scTenifoldNet Construction 🐢🐢 ...',...
%    'Callback',@callback_scTenifoldNet1);
i_addmenu(m_exp2,0,@callback_scTenifoldNet2,'scTenifoldNet Comparison 🐢🐢🐢 ...');
%uimenu(m_exp2,'Text','scTenifoldNet Comparison 🐢🐢🐢 ...',...
%    'Callback',@callback_scTenifoldNet2);

i_addmenu(m_exp2,1,@callback_scTenifoldKnk1,'scTenifoldKnk (Virtual KO) Single Gene 🐢 ...');
%uimenu(m_exp2,'Text','scTenifoldKnk (Virtual KO) Single Gene 🐢 ...',...
%    'Separator','on',...
%    'Callback',@callback_scTenifoldKnk1);

%uimenu(m_exp2,'Text','scTenifoldKnk (Virtual KO) All Genes 🐢🐢🐢 ...',...
%    'Callback',@callback_scTenifoldKnkN);

i_addmenu(m_exp,1,@callback_TrajectoryAnalysis,'Run pseudotime analysis (Monocle)...');
%uimenu(m_exp,'Text','Run pseudotime analysis (Monocle)...',...
%    'Separator','on',...
%    'Callback',@callback_TrajectoryAnalysis);

i_addmenu(m_exp,0,@DrawTrajectory,'Plot pseudotime trajectory...');
i_addmenu(m_exp,1,@callback_DetectCellularCrosstalk,'Ligand-Receptor Mediated Intercellular Crosstalk...');
i_addmenu(m_exp,0,@callback_SelectCellsByMarker,'Extract Cells by Marker(+/-) Expression...');

uimenu(m_exp,'Text','Merge Subclusters of Same Cell Type...',...
    'Callback',@MergeSubCellTypes);
i_addmenu(m_exp,0,@AnnotateSubGroup,'Annotate Cell Subgroups...');
i_addmenu(m_exp,0,@gui.callback_CellHeatMap,'Cell Heatmap...');

%uimenu(m_exp,'Text','Calculate Gene Expression Statistics...',...
%    'Callback',@callback_CalculateGeneStats);
%uimenu(m_exp,'Text','Library Size of Cell Cycle Phases...',...
%    'Callback',@callback_CellCycleLibrarySize);

i_addmenu(m_exp,1,@callback_CalculateGeneStats,'Calculate Gene Expression Statistics...');

i_addmenu(m_exp,0,@callback_CellCycleLibrarySize,'Library Size of Cell Cycle Phases...');

%uimenu(m_exp,'Text','Show HgB-genes Expression...',...
%    'Callback',@callback_ShowHgBGeneExpression);
i_addmenu(m_exp,0,@callback_ShowHgBGeneExpression,'Show HgB-genes Expression...');

%uimenu(m_exp,'Text','Show Mt-genes Expression...',...
%    'Callback',@callback_ShowMtGeneExpression);
i_addmenu(m_exp,0,@callback_ShowMtGeneExpression,'Show Mt-genes Expression...');

%uimenu(m_exp,'Text','T Cell Exhaustion Score...',...
%    'Callback',@callback_TCellExhaustionScores);
i_addmenu(m_exp,0,@callback_TCellExhaustionScores,'T Cell Exhaustion Score...');
i_addmenu(m_exp,1,@GEOAccessionToSCE,'Import Data Using GEO Accession...');
%uimenu(m_exp,'Text','Import Data Using GEO Accession...',...
%    'Separator','on',...
%    'Callback',@GEOAccessionToSCE);

%uimenu(m_exp,'Text','Merge SCEs in Workspace...',...    
%    'Callback',{@MergeSCEs,1});

i_addmenu(m_exp,0,{@MergeSCEs,1},'Merge SCEs in Workspace...');
i_addmenu(m_exp,0,{@MergeSCEs,2},'Merge SCE Data Files...');

%uimenu(m_exp,'Text','Merge SCE Data Files...',...
%    'Callback',{@MergeSCEs,2});
%uimenu(m_exp,'Text','Export Graphics as SVG File...','Callback',{@i_savefig,1});
%uimenu(m_exp,'Text','Export Graphics as Image File...','Callback',{@i_savefig,2});


i_addmenu(m_exp,0,{@i_savefig,1},'Save as SVG File...');
i_addmenu(m_exp,0,{@i_savefig,2},'Export Graphics...');

uimenu(m_exp,'Text','View Metadata...','Callback',@callback_ViewMetaData);

uimenu(m_exp,'Text','Check for Updates...',...    
    'Callback',@callback_CheckUpdates);

% handles = guihandles( FigureHandle ) ;
% guidata( FigureHandle, handles ) ;

set(FigureHandle, 'visible', 'on');
guidata(FigureHandle, sce);
setappdata(FigureHandle,'cL',cL);

set(FigureHandle,'CloseRequestFcn',@closeRequest);

if nargout > 0
    varargout{1} = FigureHandle;
end


    function i_savefig(~,~,tag)        
        if tag==1
            filter = {'*.svg'};
            [filename,filepath] = uiputfile(filter);
            if ischar(filename)
                saveas(FigureHandle,[filepath filename],'svg');
            end
        elseif tag==2
            % axx=gca;
            filter = {'*.jpg';'*.png';'*.tif';'*.pdf';'*.eps'};
            [filename,filepath] = uiputfile(filter);
            if ischar(filename)
                exportgraphics(FigureHandle,[filepath filename]);
            end
        end
    end

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

    function i_addbutton(toolbarHdl,sepTag,callbackFnc,imgFil,tooltipTxt)
        if ischar(callbackFnc) || isstring(callbackFnc)
            callbackFnc=str2func(callbackFnc);
        end
        if sepTag==1
            septag='on';
        else
            septag='off';
        end
        if toolbarHdl==1
            barhandle=UitoolbarHandle;
        else
            barhandle=defaultToolbar;
        end
        pt = uipushtool(barhandle, 'Separator', septag);
        [img, map] = imread(fullfile(mfolder, 'resources', imgFil));
        ptImage = ind2rgb(img, map);
        pt.CData = ptImage;
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
    end


% ------------------------
% Callback Functions
% ------------------------

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
    
    
    function MergeSCEs(src, ~, sourcetag)
        [requirerefresh,s]=gui.callback_MergeSCEs(src,sourcetag);
        if requirerefresh && ~isempty(s)
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c_batch_id);
            RefreshAll(src, 1, true);
            msgbox(sprintf('%s SCEs merged.', upper(s)));
        end
    end

    function AnnotateSubGroup(src, ~)
        [requirerefresh,highlightindex,newclassidenty]=gui.callback_AnnotateSubGroup(src);
        if requirerefresh && ~isempty(newclassidenty)
            disp('OK');
            sce.c(highlightindex)=sce.c(highlightindex)+newclassidenty./10;
            [c, cL] = grp2idx(sce.c);
            RefreshAll(src, 1, true);
            answer=questdlg('Save current classes to SCE.C_CLUSTER_ID?');
            if strcmp(answer,'Yes')
                sce.c_cluster_id=c;
                guidata(FigureHandle,sce);
            end
        end
    end

    function SelectCellsByQC(src, ~)
        oldn=sce.NumCells;
        oldm=sce.NumGenes;
        [requirerefresh,highlightindex]=gui.callback_SelectCellsByQC(src);
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

    function RunSeuratWorkflow(src,~)
       answer = questdlg('Run Seurat standard worflow?');
       if ~strcmp(answer, 'Yes'), return; end
       [ndim]=gui.i_choose2d3d;
       if isempty(ndim), return; end       
	   fw = gui.gui_waitbar;
       try
           [sce]=run.SeuratWorkflow(sce,ndim);
           [c, cL] = grp2idx(sce.c);
       catch
       	   gui.gui_waitbar(fw);
           return;
       end
	   gui.gui_waitbar(fw);
       RefreshAll(src, 1, true, false);
    end

    function DecontX(~,~)
        fw = gui.gui_waitbar;
        try
        [Xdecon,contamination]=run.decontX(sce);
        catch
            gui.gui_waitbar(fw);
            errordlg('Runtime error.')
            return;
        end
        gui.gui_waitbar(fw);
        figure;
        gui.i_stemscatter(sce.s,contamination);
        % zlim([0 1]);
        zlabel('Contamination rate')
        title('Ambient RNA contamination')
        answer=questdlg("Remove contamination?");
        switch answer
            case 'Yes'
                sce.X=round(Xdecon);
                guidata(FigureHandle,sce);
                helpdlg('Contamination removed.')
       end
    end

    function HarmonyPy(src, ~)
        if gui.callback_Harmonypy(src)
            sce = guidata(FigureHandle);
            [c, cL] = grp2idx(sce.c);
            RefreshAll(src, 1, true, false);
            ButtonName = questdlg('Update Saved Embedding?','', ...
                'tSNE','UMAP','PHATE','tSNE');
            methodtag=lower(ButtonName);
            if ismember(methodtag,{'tsne','umap','phate'})
                sce.struct_cell_embeddings.(methodtag)=sce.s;
            end
        end
        guidata(FigureHandle, sce);
    end


    function DoubletDetection(src, ~)        
        [isDoublet,doubletscore,methodtag,done]=gui.callback_DoubletDetection(src);
        if done && ~any(isDoublet)
            helpdlg('No doublet detected.','');
            return;
        end
        if done && any(isDoublet) && sce.NumCells==length(doubletscore)
            tmpf_doubletdetection=figure;
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
        h.Marker = para.oldMarker;
        h.SizeData = para.oldSizeData;
        colormap(para.oldColorMap);
    end

    function RenameCellType(~, ~)
        if isempty(sce.c_cell_type_tx)
            errordlg('sce.c_cell_type_tx undefined');
            return;
        end
        answer = questdlg('Rename a cell type?');
        if ~strcmp(answer, 'Yes'), return; end
        [ci, cLi] = grp2idx(sce.c_cell_type_tx);
        [indxx, tfx] = listdlg('PromptString',...
            {'Select cell type'},...
            'SelectionMode', 'single',...
            'ListString', string(cLi));
        if tfx == 1
            i = ismember(ci, indxx);
            newctype = inputdlg('New cell type', 'Rename', [1 50], cLi(ci(i)));
            if ~isempty(newctype)
                cLi(ci(i)) = newctype;
                sce.c_cell_type_tx = string(cLi(ci));
                [c, cL] = grp2idx(sce.c_cell_type_tx);
                i_labelclusters(false);
            end
        end
        guidata(FigureHandle, sce);
    end

    function EmbeddingAgain(src, ~)
        answer = questdlg('Which embedding method?', 'Select method',...
                          'tSNE', 'UMAP', 'PHATE', 'tSNE');
        if ~ismember(answer, {'tSNE', 'UMAP', 'PHATE'}), return; end
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
                '', '2000 HVGs', 'All Genes', 'Cancel', '2000 HVGs');
            switch answer2
                case 'All Genes'
                    usehvgs = false;
                case '2000 HVGs'
                    usehvgs = true;
                case {'Cancel', ''}
                    return;
            end
            [ndim]=gui.i_choose2d3d;
            if isempty(ndim), return; end
            fw = gui.gui_waitbar;
            try
                forced = true;
                sce = sce.embedcells(methodtag, forced, usehvgs, ndim);
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

    function DetermineCellTypeClusters(src, ~)
        answer = questdlg('Assign cell types to clusters automatically?',...
            '','Yes, automatically','No, manually',...
            'Cancel','Yes, automatically');
        switch answer
            case 'Yes, automatically'
                manuallyselect=false;
                bestonly=true;
            case 'No, manually'
                manuallyselect=true;
                bestonly=false;
            otherwise
                return;
        end
        speciestag = gui.i_selectspecies;
        if isempty(speciestag), return; end        
        organtag = "all";
        databasetag = "panglaodb";
        dtp = findobj(h,'Type','datatip');
        delete(dtp);
        cLdisp = cL;
        if ~manuallyselect, fw=gui.gui_waitbar_adv; end
        for i = 1:max(c)
            if ~manuallyselect
                gui.gui_waitbar_adv(fw,i/max(c));
            end
            ptsSelected = c == i;
            [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, ...
                sce.s, ptsSelected, ...
                speciestag, organtag, databasetag, bestonly);
            if isempty(Tct)
                ctxt={'Unknown'};
            else
                ctxt=Tct.C1_Cell_Type;
            end
            
            if manuallyselect
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
                [k] = dsearchn(siv, si);    % Nearest point search
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
        [h] = gui.i_gscatter3(sce.s, c, methodid, hAx);
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
            [indx, tf] = listdlg('PromptString',...
                {'Select target cluster'}, 'SelectionMode',...
                'single', 'ListString', string(c_members));
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
        [h] = gui.i_gscatter3(sce.s, c, methodid, hAx);
        title(sce.title);
        view(ax, bx);
        i_labelclusters(true);
        sce.c_cluster_id = c;
        guidata(FigureHandle, sce);
    end

    function Brush4Celltypes(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return;
        end

        answer = questdlg('Labels are not saved. Continue?');
        if ~strcmp(answer, 'Yes')
            return;
        end

        speciestag = gui.i_selectspecies;
        if isempty(speciestag), return; end  
        fw = gui.gui_waitbar;
        [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, sce.s,...
            ptsSelected, ...
            speciestag, "all", "panglaodb");
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
        
        %{
        hold on;
        if size(sce.s, 2) >= 3
            %scatter3(sce.s(ptsSelected, 1), sce.s(ptsSelected, 2),...
            %    sce.s(ptsSelected, 3), 'x');
            si = mean(sce.s(ptsSelected, :));
            text(si(:, 1), si(:, 2), si(:, 3), sprintf('%s', ctxt), ...
                'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor',...
                'w', 'EdgeColor', 'k');
        elseif size(sce.s, 2) == 2
            % scatter(sce.s(ptsSelected, 1), sce.s(ptsSelected, 2), 'x');
            si = mean(sce.s(ptsSelected, :));
            text(si(:, 1), si(:, 2), sprintf('%s', ctxt), ...
                'fontsize', 10, 'FontWeight', 'bold',...
                'BackgroundColor', 'w', 'EdgeColor', 'k');
        end
        hold off;
        %}
    
    function ShowCellStemScatter(~, ~)
        sce = guidata(FigureHandle);
        [thisc,clable]=gui.i_select1state(sce);
        if isempty(thisc), return; end        
        figure;
        gui.i_stemscatter(sce.s,grp2idx(thisc));
        zlabel(clable)
        % title('')
    end

    function ShowCellStates(src, ~)
        sce=guidata(FigureHandle);
        
        [thisc,clable,~,newpickclable]=gui.i_select1state(sce);
        
        if strcmp(clable,'Cell Cycle Phase')
            if length(unique(thisc))>1
                sce.c_cell_cycle_tx=thisc;
            end
        end
        if isempty(thisc), return; end
        if strcmp(clable,'Customized C...')
            clable=gui.i_renamec(clable,sce,newpickclable);
            sce.list_cell_attributes=[sce.list_cell_attributes,clable];
            sce.list_cell_attributes=[sce.list_cell_attributes,thisc];
        end
        [c,cL]=grp2idx(thisc);        
        sce.c=c;
        RefreshAll(src, 1, true, false);
        n=max(c);
        if n<40
            f=0.5*(n-1)./n;
            f=1+f.*(1:2:2*n);        
            cb=colorbar('Ticks',f,'TickLabels',cellstr(cL));
        else
            cb=colorbar;
        end
        cb.Label.String = clable;
        % helpdlg(clable)
        % guidata(FigureHandle, sce);
    end

    function DeleteSelectedCells(src, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No cells are selected.");
            return;
        end
        [ptsSelected,letdoit]=gui.i_expandbrushed(ptsSelected,sce);
        if ~letdoit, return; end
        
%         if ~any(ptsSelected)
%             answer = questdlg('No brushed cells. Select cells by class?');
%             if ~strcmp(answer, 'Yes'), return; end            
%             [thisc,~]=gui.i_select1class(sce);
%             if isempty(thisc), return; end
%             [ci,cLi]=grp2idx(thisc);
%             [indxx,tfx] = listdlg('PromptString',{'Select groups'},...
%                 'SelectionMode','multiple','ListString',string(cLi));
%             if tfx~=1, return; end
%             ptsSelected=ismember(ci,indxx);
%         end
        answer = questdlg('Delete cells?','','Selected', 'Unselected',...
                          'Cancel', 'Selected');
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
%        [c, cL] = grp2idx(sce.c);
%        [ax, bx] = view();
%        h = gui.i_gscatter3(sce.s, c);
%        title(sce.title);
%        view(ax, bx);

    end

    function DrawTrajectory(~, ~)
        answer = questdlg('Which method?', 'Select Algorithm', ...
            'splinefit (🐇)', 'princurve (🐢)', ...
            'splinefit (🐇)');
        if strcmp(answer, 'splinefit (🐇)')
            dim = 1;
            [t, xyz1] = i_pseudotime_by_splinefit(sce.s, dim, false);
        elseif strcmp(answer, 'princurve (🐢)')
            [t, xyz1] = i_pseudotime_by_princurve(sce.s, false);
        else
            errordlg('Invalid Option.');
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
                figure;
                i_plot_pseudotimeseries(log2(sce.X + 1), ...
                    sce.g, t, selectedg);
            case 'No'
                return
        end
        
    end

    function ClusterCellsS(src, ~)
        answer = questdlg('Cluster cells?');
        if ~strcmp(answer, 'Yes'), return; end        
        answer = questdlg('Which method?', 'Select Algorithm', ...
            'K-means 🐇', 'SnnDpc 🐢', 'K-means 🐇');
        if strcmpi(answer, 'K-means 🐇')
            methodtag = "kmeans";
        elseif strcmpi(answer, 'SnnDpc 🐢')
            methodtag = "snndpc";
        else
            return;
        end
        i_reclustercells(src, methodtag);
        guidata(FigureHandle, sce);
    end

    function ClusterCellsX(src, ~)
        answer = questdlg('Cluster cells using X?');
        if ~strcmp(answer, 'Yes')
            return
        end
        
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
                warndlg('SC3 is not working in standalone application.');
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
                    return
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
        if strcmp(state, 'off')
            dtp = findobj(h, 'Type', 'datatip');
            delete(dtp);
        else
            sce=guidata(FigureHandle);
            [thisc,~]=gui.i_select1class(sce);
            if ~isempty(thisc)
                [c,cL] = grp2idx(thisc);
                
                sce.c = c;
                RefreshAll(src, 1, true, false);                
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
            end        
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
                stxtyes = c;
            else
                [~,cLx]=grp2idx(c); 
                if isequal(cL,cLx)
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
                        return;
                end
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

%     function [para] = i_getoldsettings(src)
%         ah = findobj(src.Parent.Parent, 'type', 'Axes');
%         ha = findobj(ah.Children, 'type', 'Scatter');
%         ha1 = ha(1);
%         oldMarker = ha1.Marker;
%         oldSizeData = ha1.SizeData;
%         oldColorMap = colormap;
%         para.oldMarker = oldMarker;
%         para.oldSizeData = oldSizeData;
%         para.oldColorMap = oldColorMap;
%     end

end
