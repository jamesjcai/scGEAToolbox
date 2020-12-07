function varargout = sc_scatter_sce(sce,varargin)

import pkg.*
p = inputParser;
% validTypes = {'kmeans','kmedoids','dbscan'};
% checkType = @(x) any(validatestring(x,validTypes));
checkC = @(x) size(sce.X,2)==length(x);
addRequired(p,'sce',@(x) isa(x,'SingleCellExperiment'));
addOptional(p,'c',sce.c,checkC);
addOptional(p,'plotit',true,@islogical);
addOptional(p,'methodid',1,@isnumeric);
parse(p,sce,varargin{:});
cin=p.Results.c;
plotit=p.Results.plotit;
methodid=p.Results.plotit;

if ~isa(sce,'SingleCellExperiment')
    error('requires sce=SingleCellExperiment();');
end

[c,cL]=grp2idx(cin);

FigureHandle = figure('Name','sc_scatter_sce','visible','off');
movegui( FigureHandle, 'center' ) ; 

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
ptaddcluster.ClickedCallback = @Brushed2Cluster;

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

handles = guihandles( FigureHandle ) ; 
guidata( FigureHandle, handles ) ; 
 
set( FigureHandle, 'visible', 'on' ) ; 

if nargout > 0
    varargout{1} = FigureHandle ; 
end

% =========================
function RefreshAll(~,~)
    if ~isempty(h.ZData)
        [ax,bx]=view();
        h=i_gscatter3(sce.s,c,methodid);
        view(ax,bx);
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

function Switch2D3D(~,~)
    oldcmp=colormap();
    if isempty(h.ZData)
        h=i_gscatter3(sce.s,c,methodid);
    else
        h=i_gscatter3(sce.s(:,1:2),c,methodid);
    end
    colormap(oldcmp);
    title(sce.title)
end

function DEGene2Groups(~,~)
    answer = questdlg('Compare two batch groups (DE gene analysis)?');
    if ~strcmp(answer,'Yes'), return; end
    if isempty(sce.c_batch_id)
        warndlg("sce.c_batch_id is empty");
        return;
    end
    
    f = waitbar(0,'Please wait...');
    pause(.5); waitbar(.67,f,'Processing your data');
    T=run_mast(sce.X(:,sce.c_batch_id==1),...
            sce.X(:,sce.c_batch_id==2),sce.g);
    waitbar(1,f,'Finishing');
    pause(1); close(f);    
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
    answer = questdlg('Embedding cells?');
    if ~strcmp(answer,'Yes'), return; end
    answer = questdlg('Which method?','Select method','tSNE','Phate','UMAP','Phate');
    f = waitbar(0,'Please wait...');
    pause(.5); waitbar(.67,f,'Processing your data');
    if strcmp(answer,'tSNE')
        sce.s=sc_tsne(sce.X,3,false);
    elseif strcmp(answer,'Phate')
        sce.s=run_phate(sce.X,3,false);
    elseif strcmp(answer,'UMAP')
        sce.s=run_umap(sce.X,false);
    end
    waitbar(1,f,'Finishing');
    pause(1); close(f);
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

function Brushed2Cluster(~,~)
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
    f = waitbar(0,'Please wait...');
    pause(.5)
    waitbar(.67,f,'Processing your data');
    [Tct]=local_celltypebrushed(sce.X,sce.g,sce.s,ptsSelected,...
          speciestag,organtag,databasetag);
    ctxt=Tct.C1_Cell_Type;            
    waitbar(1,f,'Finishing');
    pause(1);
    close(f);
    
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
    answer = questdlg('Get marker genes of brushed cells?');
    if ~strcmp(answer,'Yes'), return; end  
    ptsSelected = logical(h.BrushData.');
    if ~any(ptsSelected)
        warndlg("No cells are selected.");
        return;
    end
    prompt = {'Enter number of panels (1-10)'};
    dlgtitle = 'Panel of 9 Genes';
    answer = inputdlg(prompt,dlgtitle,[1 40],{'1'});
    if isempty(answer)
        return; 
    else
        try
            numfig=str2double(answer{1});
        catch ME
            errordlg(ME.message);
            return;
        end
    end
    if ~(numfig>0 && numfig<=10)
        errordlg('Invalid number of figures');
        return;
    end
    f = waitbar(0,'Please wait...');
    pause(.5)
    waitbar(.67,f,'Processing your data');
    [markerlist]=sc_pickmarkers(sce.X,sce.g,1+ptsSelected,2);
    waitbar(1,f,'Finishing');
    pause(1);
    close(f);
    [ax,bx]=view();
    i_markergenespanel(sce.X,sce.g,sce.s,...
        markerlist,numfig,9,ax,bx,...
        'Marker Genes for Selected Cells');
    pause(2);
    export2wsdlg({'Save marker list to variable named:'},...
            {'g_markerlist'},{markerlist});
end

function ShowMarkerGene(~,~)
    answer = questdlg('Select a gene to show expression?');
    if ~strcmp(answer,'Yes'), return; end
    gsorted=sort(sce.g);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1        
        figure;
        sc_scattermarker(sce.X,sce.g,sce.s,gsorted(indx),5);
    end
end

function ShowCellStats(~,~)
    answer = questdlg('Show cell states?');
    if ~strcmp(answer,'Yes'), return; end
    
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
                if isempty(sce.c_cell_cycle_phase_tx)
                    f = waitbar(0,'Please wait...');
                    pause(.5)
                    waitbar(.67,f,'Processing your data');                
                    [cix]=run_cellcycle(sce.X,sce.g);
                    sce.c_cell_cycle_phase_tx=cix;
                    waitbar(1,f,'Finishing'); pause(1); close(f);                    
                labels = {'Save cell cycle phase to variable named:'};
                vars = {'c_cell_cycle_phase_tx'};
                values = {sce.c_cell_cycle_phase_tx};
                export2wsdlg(labels,vars,values);                
                end              
                [ci,tx]=grp2idx(sce.c_cell_cycle_phase_tx);
                ttxt=sprintf('%s|',string(tx));
            case 5 % cell type
                [ci]=grp2idx(sce.c_cell_type_tx);
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
           % colormap default
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
    if ~isempty(sce.c_cell_cycle_phase_tx)
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
            [ci,cLi]=grp2idx(sce.c_cell_cycle_phase_tx);
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
    labels = {'Save SCE to variable named:'}; 
    vars = {'sce'};
    values = {sce};
    export2wsdlg(labels,vars,values);
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
    
    f = waitbar(0,'Please wait...');
    pause(.5)
    waitbar(.67,f,'Processing your data');
    [t_mono,s_mono]=run_monocle(sce.X);
    waitbar(1,f,'Finishing');
    pause(1)
    close(f)
    

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
    f = waitbar(0,'Please wait...');
    pause(.5)
    waitbar(.67,f,'Processing your data');    
    hold off
    sce.c_cluster_id=sc_cluster_s(sce.s,k,'type',methodtag,'plotit',false);
    [c,cL]=grp2idx(sce.c_cluster_id);
    waitbar(1,f,'Finishing');
    pause(1);
    close(f);
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
    
    f = waitbar(0,'Please wait...');
    pause(.5)
    waitbar(.67,f,'Processing your data');
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
    waitbar(1,f,'Finishing');
    pause(1);
    close(f);
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
            RefreshAll;
        else
            answer = questdlg('Change current class type?');
            if strcmp(answer,'No')
                set(src,'State','off');
            elseif strcmp(answer,'Yes')                
                listitems={'c'};
                if ~isempty(sce.c_cluster_id), listitems=[listitems,'c_cluster_id']; end
                if ~isempty(sce.c_cell_type_tx), listitems=[listitems,'c_cell_type_tx']; end
                if ~isempty(sce.c_cell_cycle_phase_tx), listitems=[listitems,'c_cell_cycle_phase_tx']; end
                if ~isempty(sce.c_batch_id), listitems=[listitems,'c_batch_id']; end
                [indx,tf] = listdlg('PromptString',{'Select statistics',...
                '',''},'SelectionMode','single','ListString',listitems);
                if tf==1   
                    switch listitems{indx}
                        case 'c_cluster_id'
                            cc=sce.c_cluster_id;
                        case 'c_cell_type_tx'
                            cc=sce.c_cell_type_tx;
                        case 'c_cell_cycle_phase_tx'
                            cc=sce.c_cell_cycle_phase_tx;
                        case 'c_batch_id'
                            cc=sce.c_batch_id;
                        otherwise
                            cc=[];
                    end
                    if ~isempty(cc)                        
                        [c,cL]=grp2idx(cc);
                        set(src,'State','off');
                        RefreshAll;
                    end                    
                else
                    set(src,'State','off');
                    RefreshAll;
                    return;
                end
            else
                set(src,'State','off');
                RefreshAll;
                return;
            end
            if i_labelclusters
                set(src,'State','on');
            else                
                set(src,'State','off');
                RefreshAll;
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
    if isempty(cL)
        txt =sprintf('cluster=%d',c(idx));
    else
        txt =sprintf('cluster=%d, type=%s',c(idx),cL{c(idx)});
    end
end

function [isdone]=i_labelclusters
    isdone=false;
    if ~isempty(cL)
        answer = questdlg(sprintf('Label %d groups with index or text?',numel(cL)),...
            'Select Format','Index','Text','Text');
        if strcmp(answer,'Text')
            stxtyes=true;
        elseif strcmp(answer,'Index')
            stxtyes=false;
        else
            return;
        end
    end
    prompt = {'Enter font size (1-30):'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'15'};
    if stxtyes, definput = {'10'}; end
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    try
        fsz=round(str2double(answer{1}));
    catch
        return;
    end
    if ~(fsz>=1 && fsz<=30), return; end
    hold on
    for i=1:max(c)
        si=sce.s(c==i,:);
        si=mean(si,1);        
        if stxtyes
            stxt=sprintf('%s',cL{i});
        else
            stxt=sprintf('%d',i);
        end
        stxt=strrep(stxt,'_','\_');
        if ~isempty(h.ZData)    % size(sce.s,2)==3
            text(si(:,1),si(:,2),si(:,3),stxt,...
                'fontsize',fsz,'FontWeight','bold',...
                'BackgroundColor','w','EdgeColor','k');
        else
            text(si(:,1),si(:,2),stxt,...
                'fontsize',fsz,'FontWeight','bold',...
                'BackgroundColor','w','EdgeColor','k');
        end
    end
    hold off
    % helpdlg(sprintf('%d clusters are labelled.',numel(cL)));
    isdone=true;
end

end
