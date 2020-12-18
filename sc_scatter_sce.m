function varargout = sc_scatter_sce(sce,varargin)

if nargin<1, error('Usage: sc_scatter_sce(sce)'); end
if ~isa(sce,'SingleCellExperiment')
    error('requires sce=SingleCellExperiment();');
end

import pkg.*
p = inputParser;
% validTypes = {'kmeans','kmedoids','dbscan'};
% checkType = @(x) any(validatestring(x,validTypes));
checkCS = @(x) isempty(x) | size(sce.X,2)==length(x);
addRequired(p,'sce',@(x) isa(x,'SingleCellExperiment'));
addOptional(p,'c',sce.c,checkCS);
addOptional(p,'s',[],checkCS);
addOptional(p,'methodid',1,@isnumeric);
parse(p,sce,varargin{:});
cin=p.Results.c;
sin=p.Results.s;
methodid=p.Results.methodid;
ax=[]; bx=[];

if isempty(cin)
    sce.c=ones(size(sce.X,2),1);
else
    sce.c=cin;
end
if ~isempty(sin), sce.s=sin; end

[c,cL]=grp2idx(sce.c);

FigureHandle = figure('Name','sc_scatter_sce',...
'position',round(1.5*[0 0 560 420]),...    
'visible','off');
movegui(FigureHandle,'center');

hAx = axes('Parent',FigureHandle);
[h]=i_gscatter3(sce.s,c,methodid);
title(sce.title);

dt=datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

% defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle
defaultToolbar = findall(FigureHandle,'Type','uitoolbar');
% UitoolbarHandle2 = uitoolbar( 'Parent', FigureHandle ) ; 
% set( UitoolbarHandle2, 'Tag' , 'FigureToolBar2' , ... 
%     'HandleVisibility' , 'on' , ... 
%     'Visible' , 'on' ) ; 

UitoolbarHandle = uitoolbar( 'Parent', FigureHandle ) ; 
set( UitoolbarHandle, 'Tag' , 'FigureToolBar' , ... 
    'HandleVisibility' , 'off' , ... 
    'Visible' , 'on' ) ; 



% UitoolbarHandle = uitoolbar(FigureHandle);
pt3 = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','list.gif'));         
ptImage = ind2rgb(img,map);
pt3.CData = ptImage;
pt3.Tooltip = 'Select a gene to show expression';
pt3.ClickedCallback = @ShowMarkerGene;

pt3a = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','list2.gif'));         
ptImage = ind2rgb(img,map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Show cell states';
pt3a.ClickedCallback = @ShowCellStats;

pt3a = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-pointfig.gif'));         
ptImage = ind2rgb(img,map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Select cells by class';
pt3a.ClickedCallback = @SelectCellsByClass;

% ------------------

ptlabelclusters = uitoggletool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
             'toolbox','matlab','icons','plotpicker-scatter.gif'));
% map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
ptImage = ind2rgb(img,map);
ptlabelclusters.CData = ptImage;
ptlabelclusters.Tooltip = 'Label clusters';
ptlabelclusters.ClickedCallback = @LabelClusters;

% ------------------ clustering

ptaddcluster = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-glyplot-face.gif'));
ptImage = ind2rgb(img,map);
ptaddcluster.CData = ptImage;
ptaddcluster.Tooltip = 'Add brushed cells to a new cluster';
ptaddcluster.ClickedCallback = @Brushed2NewCluster;

ptmergecluster = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-pzmap.gif'));
ptImage = ind2rgb(img,map);
ptmergecluster.CData = ptImage;
ptmergecluster.Tooltip = 'Merge brushed cells to same cluster';
ptmergecluster.ClickedCallback = @Brushed2MergeClusters;

ptShowClu = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-geoscatter.gif'));         
ptImage = ind2rgb(img,map);
ptShowClu.CData = ptImage;
ptShowClu.Tooltip = 'Show clusters individually';
ptShowClu.ClickedCallback = @ShowClustersPop;


ptcluster = uipushtool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-dendrogram.gif'));
ptImage = ind2rgb(img,map);
ptcluster.CData = ptImage;
ptcluster.Tooltip = 'Clustering using using embedding S';
ptcluster.ClickedCallback = @ClusterCellsS;

ptcluster = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-gscatter.gif'));
ptImage = ind2rgb(img,map);
ptcluster.CData = ptImage;
ptcluster.Tooltip = 'Clustering using expression matrix X';
ptcluster.ClickedCallback = @ClusterCellsX;

% -------------

pt5 = uipushtool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','brush.gif'));
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Cell types of brushed cells';
pt5.ClickedCallback = @Brush4Celltypes;


ptclustertype = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(matlabroot,...
             'toolbox','matlab','icons','plotpicker-contour.gif'));
ptImage = ind2rgb(img,map);
ptclustertype.CData = ptImage;
ptclustertype.Tooltip = 'Cell types of clusters';
ptclustertype.ClickedCallback = @DetermineCellTypeClusters;


pt4 = uipushtool(UitoolbarHandle,'Separator','off');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-kagi.gif'));
ptImage = ind2rgb(img,map);
pt4.CData = ptImage;
pt4.Tooltip = 'Marker genes of brushed cells';
pt4.ClickedCallback = @Brush4Markers;


pt4 = uipushtool(UitoolbarHandle,'Separator','off');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-scatterhist.gif'));
ptImage = ind2rgb(img,map);
pt4.CData = ptImage;
pt4.Tooltip = 'Rename cell type';
pt4.ClickedCallback = @RenameCellType;

% --------------------------

ptpseudotime = uipushtool(defaultToolbar,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-priceandvol.gif'));
ptImage = ind2rgb(img,map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Compare Differentiation Potency';
ptpseudotime.ClickedCallback = @ComparePotency;

ptpseudotime = uipushtool(defaultToolbar,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-arxtimeseries.gif'));
ptImage = ind2rgb(img,map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Run pseudotime analysis (Monocle)';
ptpseudotime.ClickedCallback = @RunTrajectoryAnalysis;

ptpseudotime = uipushtool(defaultToolbar,...
    'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-comet.gif'));
ptImage = ind2rgb(img,map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Plot pseudotime trajectory';
ptpseudotime.ClickedCallback = @DrawTrajectory;


pt4 = uipushtool(defaultToolbar,'Separator','on');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-boxplot.gif'));
ptImage = ind2rgb(img,map);
pt4.CData = ptImage;
pt4.Tooltip = 'Compare 2 groups (DE analysis)';
pt4.ClickedCallback = @DEGene2Groups;

ptpseudotime = uipushtool(defaultToolbar,...
    'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-andrewsplot.gif'));
ptImage = ind2rgb(img,map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Function enrichment of HVG genes';
ptpseudotime.ClickedCallback = {@callback_GSEA_HVGs,sce.X,sce.g};

pt2 = uipushtool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-qqplot.gif'));  
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delete selected cells';
pt2.ClickedCallback = @DeleteSelectedCells;

pt = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','export.gif'));         
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Export & save data';
pt.ClickedCallback = @SaveX;


pt5 = uipushtool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-geobubble.gif'));
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Embedding';
pt5.ClickedCallback = @EmbeddingAgain;


pt5 = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-image.gif'));  % plotpicker-pie
%map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Switch 2D/3D';
pt5.ClickedCallback = @Switch2D3D;



pt5pickcolr = uipushtool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-compass.gif'));  % plotpicker-pie
%map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
ptImage = ind2rgb(img,map);
pt5pickcolr.CData = ptImage;
pt5pickcolr.Tooltip = 'Switch color maps';
pt5pickcolr.ClickedCallback = @callback_PickColormap;

pt5 = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-geobubble2.gif'));
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Refresh';
pt5.ClickedCallback = @RefreshAll;


add_3dcamera(defaultToolbar,'AllCells');

% handles = guihandles( FigureHandle ) ; 
% guidata( FigureHandle, handles ) ; 
 
set( FigureHandle, 'visible', 'on' ) ; 

if nargout > 0
    varargout{1} = FigureHandle ; 
end



% =========================
function RefreshAll(~,~)
    if size(sce.s,2)>2
        if ~isempty(h.ZData)
            %[ax,bx]=view();
            h=i_gscatter3(sce.s,c,methodid);
            %view(ax,bx);
        else
            h=i_gscatter3(sce.s(:,1:2),c,methodid);
        end       
    else
        h=i_gscatter3(sce.s(:,1:2),c,methodid);        
    end
    title(sce.title)
    
%    h=i_gscatter3(sce.s,c,methodid);
%    title(sce.title)
    
    ptlabelclusters.State='off';
    UitoolbarHandle.Visible='off';
    UitoolbarHandle.Visible='on';
    
    %legend off
    %colorbar off    
end

function ComparePotency(~,~)
    if ~isempty(sce.list_cell_attributes{2}) && ~isempty(sce.c_cell_type_tx)
        x=sce.list_cell_attributes{2};
        y=sce.c_cell_type_tx;
    figure;
    pkg.i_violinplot_groupordered(x,y);
    if ~isempty(sce.c_batch_id)
        figure;
        h1=subplot(1,2,1);
        i=sce.c_batch_id==1;
        pkg.i_violinplot_groupordered(x(i),y(i));
        [aabb]=ylim();
        h2=subplot(1,2,2);
        i=sce.c_batch_id==2;
        pkg.i_violinplot_groupordered(x(i),y(i));
        [ccdd]=ylim();
        ylim(h1,[min(aabb(1),ccdd(1)) max(aabb(2),ccdd(2))]);
        ylim(h2,[min(aabb(1),ccdd(1)) max(aabb(2),ccdd(2))]);
    end
    end
end

function Switch2D3D(~,~)
    %oldcmp=colormap();
    if isempty(h.ZData)   % current 2 D
        if ~(size(sce.s,2)>2)
            helpdlg('Canno swith to 3-D. SCE.S is 2-D');
            return;
        end        
        h=i_gscatter3(sce.s,c,methodid);
        if ~isempty(ax) && ~isempty(bx) && ~any([ax bx]==0)
            view(ax,bx);
        else
            view(3);
        end        
    else                 % current 3D do following
        [ax,bx]=view();
        answer = questdlg('Which view to be used to project cells?','',...
            'Current View','Default View','Cancel','Current View');
        if strcmp(answer,'Cancel'), return; end
        if strcmp(answer,'Default View')
            h=i_gscatter3(sce.s(:,1:2),c,methodid);
        else
            sx=pkg.i_3d2d(sce.s,ax,bx);
            h=i_gscatter3(sx(:,1:2),c,methodid);
        end    
    end
    title(sce.title)
end

function RenameCellType(~,~)
    if isempty(sce.c_cell_type_tx)
        errordlg('sce.c_cell_type_tx undefined');
        return;
    end
    answer = questdlg('Rename a cell type?');
    if ~strcmp(answer,'Yes'), return; end
    [ci,cLi]=grp2idx(sce.c_cell_type_tx);
    [indxx,tfx] = listdlg('PromptString',{'Select a cell type',...
    '',''},'SelectionMode','single','ListString',string(cLi));
    if tfx==1
        i=ismember(ci,indxx);
        newctype=inputdlg('New cell type','Rename',[1 50],cLi(ci(i)));
        if ~isempty(newctype)
            cLi(ci(i))=newctype;
            sce.c_cell_type_tx=cLi(ci);
            [c,cL]=grp2idx(sce.c_cell_type_tx);
            i_labelclusters(false);
        end
    end
    
end

function DEGene2Groups(~,~)
    answer = questdlg('Compare two batch groups (DE gene analysis)?');
    if ~strcmp(answer,'Yes'), return; end
    if isempty(sce.c_batch_id)
        warndlg("sce.c_batch_id is empty");
        return;
    end
    answer = questdlg('Which method?','Select Method','Wilcoxon rank-sum test','MAST','Wilcoxon rank-sum test');
    if strcmpi(answer,'Wilcoxon rank-sum test')
        methodtag="ranksum";
    elseif strcmpi(answer,'MAST')
        methodtag="mast";
    else
        return;
    end
    fw=pkg.gui_waitbar;
    switch methodtag
        case 'ranksum'
            T=sc_deg(sce.X(:,sce.c_batch_id==1),...
                    sce.X(:,sce.c_batch_id==2),sce.g);
        case 'mast'
            T=run_mast(sce.X(:,sce.c_batch_id==1),...
                    sce.X(:,sce.c_batch_id==2),sce.g);
    end
    pkg.gui_waitbar(fw);
    labels = {'Save DE results T to variable named:'}; 
    vars = {'T'}; values = {T};
    msgfig=export2wsdlg(labels,vars,values);
    uiwait(msgfig);
    answer = questdlg('Violin plots?');
    if strcmp(answer,'Yes')
        figure;
        for k=1:16
            subplot(4,4,k)
            i=sce.g==T.gene(k);
            pkg.i_violinplot(log2(1+sce.X(i,:)),...
                sce.c_batch_id);
            title(T.gene(k));
            ylabel('log2(UMI+1)')
        end
    end
end

function EmbeddingAgain(~,~)
    %answer = questdlg('Embedding cells?');
    %if ~strcmp(answer,'Yes'), return; end
    answer = questdlg('Which method?','Select method','tSNE','UMAP','PHATE','tSNE');
    fw=pkg.gui_waitbar;
    if strcmp(answer,'tSNE')
        sce.s=sc_tsne(sce.X,3,false);
    elseif strcmp(answer,'UMAP')
        sce.s=run_umap(sce.X,2,false,false);
    elseif strcmp(answer,'PHATE')
        sce.s=run_phate(sce.X,3,false);       
    end
    pkg.gui_waitbar(fw);
    RefreshAll;
end

function DetermineCellTypeClusters(~,~)
    answer = questdlg('Label cell type of clusters?');
    if ~strcmp(answer,'Yes'), return; end    
    answer = questdlg('Which species?','Select Species','Mouse','Human','Mouse');

    if ~strcmp(answer,'Human')
        speciestag="human";
    else
        speciestag="mouse";
    end
    organtag="all";
    
    answer = questdlg('Which marker database?','Select Database','PanglaoDB','clustermole','PanglaoDB');
    if strcmpi(answer,'clustermole')
        databasetag="clustermole";
    else
        databasetag="panglaodb";
    end    
    
for i=1:max(c) 
    ptsSelected=c==i;
    [Tct]=local_celltypebrushed(sce.X,sce.g,sce.s,ptsSelected,...
          speciestag,organtag,databasetag);
    ctxt=Tct.C1_Cell_Type;
    
    [indx,tf] = listdlg('PromptString',{'Select cell type',...
    '',''},'SelectionMode','single','ListString',ctxt);
    if tf==1
        ctxt=Tct.C1_Cell_Type{indx};
    else
        return;
    end
    hold on
    cL{i}=ctxt;
    ctxt=strrep(ctxt,'_','\_');            
    if size(sce.s,2)>=3
            si=mean(sce.s(ptsSelected,:));
            text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
                 'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
    elseif size(sce.s,2)==2
            si=mean(sce.s(ptsSelected,:));
            text(si(:,1),si(:,2),sprintf('%s',ctxt),...
                 'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
    end
    hold off
end
sce.c_cell_type_tx=string(cL(c));
export2wsdlg({'Save cell type list to variable named:'},...
    {'c_celltype'},{sce.c_cell_type_tx});
end

function Brushed2NewCluster(~,~)
    answer = questdlg('Make a new cluster out of brushed cells?');
    if ~strcmp(answer,'Yes'), return; end  
    ptsSelected = logical(h.BrushData.');
    if ~any(ptsSelected)
        warndlg("No cells are selected.");
        return;
    end
    c(ptsSelected)=max(c)+1;
    [c,cL]=grp2idx(c);
    sce.c=c;
    [ax,bx]=view();
    [h]=i_gscatter3(sce.s,c,methodid);
    title(sce.title)
    view(ax,bx);
    i_labelclusters(true);
    answer = questdlg('Update sce.c_cluster_id?');
    if strcmp(answer,'Yes')
        sce.c_cluster_id=c;
    else
        return;
    end     
end

function Brushed2MergeClusters(~,~)
    answer = questdlg('Merge brushed cells into one cluster?');
    if ~strcmp(answer,'Yes'), return; end  
    ptsSelected = logical(h.BrushData.');
    if ~any(ptsSelected)
        warndlg("No cells are brushed");
        return;
    end
    c_members=unique(c(ptsSelected));
    if numel(c_members)==1
        warndlg("All brushed cells are in one cluster");
        return;
    else
        [indx,tf] = listdlg('PromptString',{'Select target cluster',...
        '',''},'SelectionMode','single','ListString',string(c_members));
        if tf==1
            c_target=c_members(indx);
        else
            return;
        end
    end    
    c(ismember(c,c_members))=c_target;
    [c,cL]=grp2idx(c);
    sce.c=c;
    [ax,bx]=view();
    [h]=i_gscatter3(sce.s,c,methodid);
    title(sce.title)
    view(ax,bx);
    i_labelclusters(true);
    answer = questdlg('Update sce.c_cluster_id?');
    if strcmp(answer,'Yes')
        sce.c_cluster_id=c;
    else
        return;
    end     
end


function Brush4Celltypes(~,~)
    answer = questdlg('Label cell type of brushed cells?');
    if ~strcmp(answer,'Yes'), return; end
    ptsSelected = logical(h.BrushData.');
    if ~any(ptsSelected)
        warndlg("No cells are selected.");
        return;
    end
    
    answer = questdlg('Which species?','Select Species','Mouse','Human','Mouse');
    if ~strcmp(answer,'Human')
        speciestag="human";
    else
        speciestag="mouse";
    end
    organtag="all";
    
    answer = questdlg('Which marker database?','Select Database','PanglaoDB','clustermole','PanglaoDB');
    if strcmpi(answer,'clustermole')
        databasetag="clustermole";
    else
        databasetag="panglaodb";
    end    
    fw=pkg.gui_waitbar;
    [Tct]=local_celltypebrushed(sce.X,sce.g,sce.s,ptsSelected,...
          speciestag,organtag,databasetag);
    ctxt=Tct.C1_Cell_Type;            
    pkg.gui_waitbar(fw);
   
        [indx,tf] = listdlg('PromptString',{'Select cell type',...
        '',''},'SelectionMode','single','ListString',ctxt);
        if tf==1 
            ctxt=Tct.C1_Cell_Type{indx};
        else
            return;
        end
            
            hold on
            ctxt=strrep(ctxt,'_','\_');            
            if size(sce.s,2)>=3
                    scatter3(sce.s(ptsSelected,1),sce.s(ptsSelected,2),sce.s(ptsSelected,3),'x');
                    si=mean(sce.s(ptsSelected,:));
                    text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
                         'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            elseif size(sce.s,2)==2
                    scatter(sce.s(ptsSelected,1),sce.s(ptsSelected,2),'x')                    
                    si=mean(sce.s(ptsSelected,:));
                    text(si(:,1),si(:,2),sprintf('%s',ctxt),...
                         'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            end
            hold off
end

function Brush4Markers(~,~)
    [ax,bx]=view();
    ptsSelected = logical(h.BrushData.');
    if ~any(ptsSelected)
        warndlg("No cells are selected.");
        return;
    end
    if isscalar(unique(c))
        methodtag=1;
    else
        answer = questdlg('Select brushed cell group?');
        if strcmp(answer,'Yes')
            if isscalar(unique(c(ptsSelected)))
                methodtag=2;
            else
                errordlg('More than one group of brushed cells');
                return;
            end
        elseif strcmp(answer,'No')
            methodtag=1;
        else
            return;
        end
    end
    fw=pkg.gui_waitbar;
    switch methodtag
        case 1
            [markerlist,A]=sc_pickmarkers(sce.X,sce.g,1+ptsSelected,2);
        case 2
            [markerlist,A]=sc_pickmarkers(sce.X,sce.g,c,unique(c(ptsSelected)));
    end
    pkg.gui_waitbar(fw);
    % assignin('base','A',A);
    [numfig]=pkg.gui_inputdlg;
    fw=pkg.gui_waitbar;
    htmlfilename=cL{unique(c(ptsSelected))};
    pkg.i_markergeneshtml(sce,markerlist,numfig,[ax bx],htmlfilename);
    pkg.gui_waitbar(fw);
%     pause(2);
%     export2wsdlg({'Save marker list to variable named:'},...
%             {'g_markerlist'},{markerlist});
end

function ShowMarkerGene(~,~)
    [axx,bxx]=view();
    if any([axx,bxx]==0), axx=ax; bxx=bx; end
%     answer = questdlg('Select a gene to show expression?');
%     if ~strcmp(answer,'Yes'), return; end
    gsorted=sort(sce.g);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1     
        figure;
        [h1]=sc_scattermarker(sce.X,sce.g,sce.s,gsorted(indx),5);
        view(h1,axx,bxx);
    end
end

function ShowCellStats(~,~)
%     answer = questdlg('Show cell states?');
%     if ~strcmp(answer,'Yes'), return; end
    
    listitems={'Library Size','Mt-reads Ratio',...
        'Mt-genes Expression','Cell Cycle Phase',...
        'Cell Type','Cluster ID','Batch ID'};
    for k=1:2:length(sce.list_cell_attributes)
        listitems=[listitems,sce.list_cell_attributes{k}];
    end
    [indx,tf] = listdlg('PromptString',{'Select statistics',...
    '',''},'SelectionMode','single','ListString',listitems);
    if tf==1        
        switch indx
            case 1
                ci=sum(sce.X);
                ttxt="Library Size";
            case 2
                i=startsWith(sce.g,'mt-','IgnoreCase',true);
                lbsz=sum(sce.X,1);
                lbsz_mt=sum(sce.X(i,:),1);
                ci=lbsz_mt./lbsz;
                ttxt="mtDNA%";
            case 3
                idx=startsWith(sce.g,'mt-','IgnoreCase',true);
                n=sum(idx);
                if n>0
                    [ax,bx]=view();
                    if n<=9
                        i_markergenespanel(sce.X,sce.g,sce.s,...
                            sce.g(idx),[],9,ax,bx,'Mt-genes');
                    else
                        i_markergenespanel(sce.X,sce.g,sce.s,...
                            sce.g(idx),[],16,ax,bx,'Mt-genes');
                    end
                else
                    warndlg('No mt-genes found');
                end
                return;
            case 4   % "Cell Cycle Phase";
                if isempty(sce.c_cell_cycle_tx)
                    fw=pkg.gui_waitbar;               
                    [cix]=run_cellcycle(sce.X,sce.g);
                    sce.c_cell_cycle_tx=cix;
                    pkg.gui_waitbar(fw);                    
                end                
                [ci,tx]=grp2idx(sce.c_cell_cycle_tx);
                ttxt=sprintf('%s|',string(tx));
            case 5 % cell type
                ci=sce.c_cell_type_tx;
            case 6 % cluster id              
                ci=sce.c_cluster_id;
            case 7 % batch id
                ci=sce.c_batch_id;
            otherwise % other properties
                ttxt=sce.list_cell_attributes{indx-7};
                ci=sce.list_cell_attributes{indx-7+1};
        end
        if isempty(ci)
            errordlg("Undefined classification");
            return;
        end
        sces=sce.s;
        if isempty(h.ZData)
            sces=sce.s(:,1:2);
        end
            [ax,bx]=view();     
            h=i_gscatter3(sces,ci,1);            
            view(ax,bx);
            title(sce.title);
            if indx==4
                hc=colorbar;
                hc.Label.String=ttxt;
            else
                colorbar off
            end
    answer = questdlg('Update sce.c?');
    if strcmp(answer,'Yes')
       %[c,cL]=grp2idx(ci);
       sce.c=ci;
    end
    end
end

function SelectCellsByClass(~,~)
    answer = questdlg('Select cells by class?');
    if ~strcmp(answer,'Yes'), return; end
    
    listitems={'Custom input (C)'};
    if ~isempty(sce.c_cluster_id)
        listitems=[listitems,'Cluster ID'];
    end
    if ~isempty(sce.c_cell_type_tx)
        listitems=[listitems,'Cell Type'];
    end
    if ~isempty(sce.c_cell_cycle_tx)
        listitems=[listitems,'Cell Cycle Phase'];
    end
    if ~isempty(sce.c_batch_id)
        listitems=[listitems,'Batch ID'];
    end
    
    [indx,tf] = listdlg('PromptString',{'Select statistics','',''},...    
    'SelectionMode','single','ListString',listitems);
    if tf~=1, return; end
    switch listitems{indx}
        case 'Custom input (C)'
            ci=c; cLi=cL;
        case 'Batch ID'
            [ci,cLi]=grp2idx(sce.c_batch_id);
        case 'Cluster ID'
            [ci,cLi]=grp2idx(sce.c_cluster_id);
        case 'Cell Type'
            [ci,cLi]=grp2idx(sce.c_cell_type_tx);
        case 'Cell Cycle Phase'
            [ci,cLi]=grp2idx(sce.c_cell_cycle_tx);
    end
    [indxx,tfx] = listdlg('PromptString',{'Select groups',...
    '',''},'SelectionMode','multiple','ListString',string(cLi));
    if tfx==1
        i=ismember(ci,indxx);
        [ax,bx]=view();
        scex=selectcells(sce,i);
        scex.c=cLi(ci(i));
        sc_scatter_sce(scex);
        view(ax,bx);
    end
end

function DeleteSelectedCells(~,~)
    ptsSelected = logical(h.BrushData.');    
    if ~any(ptsSelected)
        warndlg("No cells are selected.");
        return;
    end    
    answer = questdlg('Delete cells?','',...
       'Selected','Unselected','Cancel','Selected');
    if strcmp(answer,'Cancel'), return; end 
    if strcmp(answer,'Unselected')
        ptsSelected=~ptsSelected; 
    end
    answer2 = questdlg(sprintf('Delete %s cells?',...
        lower(answer)));
    if ~strcmp(answer2,'Yes'), return; end
    sce=rmcells(sce,ptsSelected);
    [c,cL]=grp2idx(sce.c);
    [ax,bx]=view();
    h=i_gscatter3(sce.s,c);
    title(sce.title);
    view(ax,bx);
end

function SaveX(~,~)
    answer = questdlg('Export & save data?');
    if ~strcmp(answer,'Yes'), return; end     
    labels = {'Save SCE to variable named:',...
        'Save SCE.X to variable named:',...
        'Save SCE.g to variable named:',...
        'Save SCE.S to variable named:'}; 
    vars = {'sce','X','genelist','s'};
    values = {sce,sce.X,sce.g,sce.s};
    export2wsdlg(labels,vars,values,...
        'Save Data to Workspace',...
        logical([1 0 0 0]),{@smhelp});
         
    function smhelp
        helpdlg({'Select one or both check boxes.',...
                 'Change the variable names, if desired,',...
                 'and then click OK.'});
    end
end

function DrawTrajectory(~,~)
        answer = questdlg('Which method?','Select Algorithm',...
            'splinefit (fast)','princurve (slow)',...
            'splinefit (fast)');
    if strcmp(answer,'splinefit (fast)')
        dim=1;
        [t,xyz1]=i_pseudotime_by_splinefit(sce.s,dim,false);    
    else
        [t,xyz1]=i_pseudotime_by_princurve(sce.s,false);
    end
        hold on
        if size(xyz1,2)>=3
            plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'-r','linewidth',2);
            text(xyz1(1,1),xyz1(1,2),xyz1(1,3),'Start',...
              'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            text(xyz1(end,1),xyz1(end,2),xyz1(end,3),'End',...
              'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');      
        elseif size(xyz1,2)==2
            plot(xyz1(:,1),xyz1(:,2),'-r','linewidth',2);
            text(xyz1(1,1),xyz1(1,2),'Start',...
              'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            text(xyz1(end,1),xyz1(end,2),'End',...
              'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');        
        end
        hold off

        labels = {'Save expression X to variable named:',...
                  'Save pseudotime T to variable named:',...
                  'Save embedding S to variable named:'}; 
        vars = {'X_psexplorer','t_psexplorer','s_psexplorer'};
        values = {sce.X, t, sce.s};
        msgfig=export2wsdlg(labels,vars,values);
        %         assignin('base',sprintf('psexplorerT%d',...
        %                  psexplorer_timeid),t);
        uiwait(msgfig)
        answer = questdlg('View expression of selected genes', ...
            'Pseudotime Function', ...
            'Yes','No','Yes');
        switch answer
            case 'Yes'
                r=corr(t,sce.X.','type','spearman'); % Calculate linear correlation between gene expression profile and T
                [~,idxp]= maxk(r,4);  % Select top 4 positively correlated genes
                [~,idxn]= mink(r,3);  % Select top 3 negatively correlated genes
                selectedg=sce.g([idxp idxn]);        
                figure;
                i_plot_pseudotimeseries(log2(sce.X+1),...
                    sce.g,t,selectedg);
            case 'No'
                return;
        end
        
end

function RunTrajectoryAnalysis(~,~)
    answer = questdlg('Run pseudotime analysis (Monocle)?');
    if ~strcmp(answer,'Yes'), return; end
    
    fw=pkg.gui_waitbar;
    [t_mono,s_mono]=run_monocle(sce.X);
    pkg.gui_waitbar(fw);

        answer = questdlg('View Monocle DDRTree?', ...
            'Pseudotime View', ...
            'Yes','No','Yes');
        switch answer
            case 'Yes'
            [ax,bx]=view();
            cla(hAx);
            sce.s=s_mono; sce.c=t_mono;
            [c,cL]=grp2idx(sce.c);
            h=i_gscatter3(sce.s,c);
            title(sce.title)
            view(ax,bx);
            hc=colorbar;
            hc.Label.String='Pseudotime';            
        end

    labels = {'Save pseudotime T to variable named:',...
              'Save embedding S to variable named:'}; 
    vars = {'t_mono','s_mono'};
    values = {t_mono, s_mono};
    export2wsdlg(labels,vars,values);          
        
end

function ClusterCellsS(~,~)
    answer = questdlg('Cluster cells?');
    if ~strcmp(answer,'Yes'), return; end

    answer = questdlg('Which method?','Select Algorithm',...
        'kmeans','snndpc','kmeans');
    if strcmpi(answer,'kmeans')
        methodtag="kmeans";
    else
        methodtag="snndpc";
    end 
    
    prompt = {'Enter number of cluster (K):'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'10'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    try
        k=str2double(answer{1});
    catch
        return;
    end    
    fw=pkg.gui_waitbar;  
    hold off
    sce.c_cluster_id=sc_cluster_s(sce.s,k,'type',methodtag,'plotit',false);
    [c,cL]=grp2idx(sce.c_cluster_id);
    pkg.gui_waitbar(fw);
    [ax,bx]=view();
    sce.c=c;
    h=i_gscatter3(sce.s,c);
    view(ax,bx);
    title(sce.title)
    answer = questdlg('Label clusters?');
    if strcmp(answer,'Yes')
        i_labelclusters;
    end
end

function ClusterCellsX(~,~)
    answer = questdlg('Cluster cells using X?');
    if ~strcmp(answer,'Yes'), return; end
    methodtagv={'simlr','soptsc','sc3','sinnlrr'};
    [indx,tf] = listdlg('PromptString',{'Select clustering program',...
    '',''},'SelectionMode','single',...
    'ListString',methodtagv);
    if tf==1
        prompt = {'Enter number of cluster (K):'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'10'};
        answerk = inputdlg(prompt,dlgtitle,dims,definput);
        try
            k=str2double(answerk{1});
        catch
            return;
        end
        methodtag=methodtagv{indx};
    else
        retrun;
    end
    
    fw=pkg.gui_waitbar;
    try
        [sce.c_cluster_id]=sc_cluster_x(sce.X,k,'type',methodtag);
    catch ME
        close(f);        
        errordlg(sprintf('%s: %s',...
            ME.identifier,ME.message));
        return;
    end
    [c,cL]=grp2idx(sce.c_cluster_id);
    sce.c=c;
    pkg.gui_waitbar(fw);
    % hold off
    delete(h);
    h=i_gscatter3(sce.s,c);
    title(sce.title)
    
    labels = {'Save clusterid C to variable named:'}; 
    vars = {sprintf('c_clusterid_%s',methodtag)};
    values = {sce.c_cluster_id};
    msgfig=export2wsdlg(labels,vars,values);
    uiwait(msgfig);
    answer = questdlg('Label clusters?');
    if strcmp(answer,'Yes')
        i_labelclusters;
    end
end

function LabelClusters(src,~)
        state = src.State;
        if strcmp(state,'off')
            dtp = findobj(h,'Type','datatip');
            delete(dtp);
        else
            answer = questdlg('Change current class type?');
            if strcmp(answer,'No')
                set(src,'State','off');
            elseif strcmp(answer,'Yes')                
                listitems={'c'};
                if ~isempty(sce.c_cluster_id), listitems=[listitems,'c_cluster_id']; end
                if ~isempty(sce.c_cell_type_tx), listitems=[listitems,'c_cell_type_tx']; end
                if ~isempty(sce.c_cell_cycle_tx), listitems=[listitems,'c_cell_cycle_tx']; end
                if ~isempty(sce.c_batch_id), listitems=[listitems,'c_batch_id']; end
                [indx,tf] = listdlg('PromptString',{'Select statistics',...
                '',''},'SelectionMode','single','ListString',listitems);
                if tf==1   
                    switch listitems{indx}
                        case 'c_cluster_id'
                            cc=sce.c_cluster_id;
                        case 'c_cell_type_tx'
                            cc=sce.c_cell_type_tx;
                        case 'c_cell_cycle_tx'
                            cc=sce.c_cell_cycle_tx;
                        case 'c_batch_id'
                            cc=sce.c_batch_id;
                        otherwise
                            cc=[];
                    end
                    if ~isempty(cc)                        
                        [c,cL]=grp2idx(cc);
                        set(src,'State','off');
                    end                    
                else
                    set(src,'State','off');
                    return;
                end
            else
                set(src,'State','off');
                return;
            end
            if i_labelclusters
                set(src,'State','on');
            else                
                set(src,'State','off');
            end
        end
end

function ShowClustersPop(~,~)
    answer = questdlg('Show clusters in new figures?');
    if ~strcmp(answer,'Yes'), return; end
    
    cmv=1:max(c);
    idxx=cmv;
    [cmx]=countmember(cmv,c);
    answer = questdlg('Sort by size of cell groups?');
    if strcmpi(answer,'Yes')        
        [~,idxx]=sort(cmx,'descend');
    end
    sces=sce.s;
    if isempty(h.ZData), sces=sce.s(:,1:2); end
    figure;
    for k=1:9
        if k<=max(c)
            subplot(3,3,k);
            i_gscatter3(sces,c,3,cmv(idxx(k)));
            title(sprintf('%s\n[%d cells (%.2f%%)]',...
                cL{idxx(k)},cmx(idxx(k)),100*cmx(idxx(k))/length(c)));
        end
    end
    
    if ceil(max(c)/9)==2
        figure;
        for k=1:9
            kk=k+9;
            if kk<=max(c)
                subplot(3,3,k);
                i_gscatter3(sces,c,3,cmv(idxx(kk)));
                title(sprintf('%s\n[%d cells (%.2f%%)]',...
                    cL{idxx(kk)},cmx(idxx(kk)),100*cmx(idxx(kk))/length(c)));
            end
        end
    end
    if ceil(max(c)/9)>2
        warndlg('Group(s) #18 and above are not displayed');
    end
end

function [txt]=i_myupdatefcnx(~,event_obj)
    % pos = event_obj.Position;
    idx = event_obj.DataIndex;
    txt=cL(c(idx));
end

function [isdone]=i_labelclusters(notasking)
    if nargin<1, notasking=false; end
    isdone=false;
    if ~isempty(cL)
        if notasking
            stxtyes=c;
        else
            answer = questdlg(sprintf('Label %d groups with index or text?',numel(cL)),...
                'Select Format','Index','Text','Text');
            if strcmp(answer,'Text')
                stxtyes=cL(c);
            elseif strcmp(answer,'Index')
                stxtyes=c;
            else
                return;
            end            
        end
        
        dtp = findobj(h,'Type','datatip');
        delete(dtp);
            
        row = dataTipTextRow('',stxtyes);
        h.DataTipTemplate.DataTipRows = row;
        for i=1:max(c)
            idx=find(c==i);
            siv=sce.s(idx,:);
            si=mean(siv,1);
            [k]=dsearchn(siv,si);
            dtp=datatip(h,'DataIndex',idx(k));            
        end
        isdone=true;
    end
end
end