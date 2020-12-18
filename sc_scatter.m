function sc_scatter(X,genelist,s,c,methodid)

if isa(X,'SingleCellExperiment')
    % matlab.lang.OnOffSwitchState.on
    sce=X;
    X=sce.X;
    genelist=sce.genelist;
    s=sce.s;
    c=sce.c;
end
if nargin<5, methodid=1; end
if nargin<4 || isempty(c), c=ones(size(X,2),1); end
if nargin<3 || isempty(s), s=randn(size(X,2),3); end
if nargin<2 || isempty(genelist), genelist=string((1:size(X,1))'); end

sce=SingleCellExperiment(X,genelist,s,c);
sc_scatter_sce(sce);

return;

if length(genelist)~=size(X,1)
    error('ERROR: length(genelist)!=size(X,1)')
end
if size(s,1)~=size(X,2)
    error('ERROR: size(s,1)!=size(X,2)')
end
if size(c,1)~=size(X,2)
    error('ERROR: size(c,1)!=size(X,2)')
end

c_cell_idx=(1:size(X,2))';
c_cell_cycle_phase=[];
c_cell_type=[];
c_cluster_id=[];
[c,cL]=grp2idx(c);

%c_custom_id=string(cL(c));
%kc=numel(unique(c));
 
% generate a new figure . 
FigureHandle = figure('Name','sc_scatter','visible','off');
movegui( FigureHandle, 'center' ) ; 


hAx = axes('Parent',FigureHandle);
[h]=i_gscatter3(s,c,methodid);
title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)));

dt=datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

% generate a uitoolbar . 
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

ptShowClu = uipushtool(UitoolbarHandle,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-geoscatter.gif'));         
ptImage = ind2rgb(img,map);
ptShowClu.CData = ptImage;
ptShowClu.Tooltip = 'Show clusters individually';
ptShowClu.ClickedCallback = @ShowClustersPop;

% ------------------ clustering

ptaddcluster = uipushtool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-pzmap.gif'));
ptImage = ind2rgb(img,map);
ptaddcluster.CData = ptImage;
ptaddcluster.Tooltip = 'Add brushed cells to a new cluster';
ptaddcluster.ClickedCallback = @Brushed2Cluster;

ptcluster = uipushtool(UitoolbarHandle,'Separator','off');
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


pt4 = uipushtool(UitoolbarHandle,'Separator','on');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-kagi.gif'));
ptImage = ind2rgb(img,map);
pt4.CData = ptImage;
pt4.Tooltip = 'Marker genes of brushed cells';
pt4.ClickedCallback = @Brush4Markers;


ptpseudotime = uipushtool(UitoolbarHandle,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-arxtimeseries.gif'));
ptImage = ind2rgb(img,map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Run pseudotime analysis (Monocle)';
ptpseudotime.ClickedCallback = @RunTrajectoryAnalysis;

ptpseudotime = uipushtool(UitoolbarHandle,...
    'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'resources','plotpicker-comet.gif'));
ptImage = ind2rgb(img,map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Plot pseudotime trajectory';
ptpseudotime.ClickedCallback = @DrawTrajectory;


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
            'resources','plotpicker-geobubble2.gif'));
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Refresh';
pt5.ClickedCallback = @RefreshAll;

add_3dcamera(UitoolbarHandle,'AllCells');

%warning off
%WinOnTop(hFig,true);
%warning on

handles = guihandles( FigureHandle ) ; 
guidata( FigureHandle, handles ) ; 
 
set( FigureHandle, 'visible', 'on' ) ; 

if nargout > 0
    varargout{1} = FigureHandle ; 
end

% =========================
function RefreshAll(~,~)
    %cla(hAx);    
    delete(h);
    % pause(.5);
    h=i_gscatter3(s,c,methodid);
    title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
    ptlabelclusters.State='off';
    UitoolbarHandle.Visible='off';
    UitoolbarHandle.Visible='on';
    legend off
    colorbar off
end

function EmbeddingAgain(~,~)
    answer = questdlg('Embedding cells?');
    if ~strcmp(answer,'Yes'), return; end
    answer = questdlg('Which method?','Select method','tSNE','Phate','UMAP','Phate');
         f = waitbar(0,'Please wait...');
     pause(.5)
     waitbar(.67,f,'Processing your data');
    if strcmp(answer,'tSNE')
        s=sc_tsne(X,3,false);
    elseif strcmp(answer,'Phate')
        s=run_phate(X,3,false);
    elseif strcmp(answer,'UMAP')
        s=run_umap(X,false);
    end
        waitbar(1,f,'Finishing');
    pause(1);
    close(f);
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
    
%     f = waitbar(0,'Please wait...');
%     pause(.5)
%     waitbar(.67,f,'Processing your data');
    
for i=1:max(c)    
%     Xi=X(:,c==i);    
%     [Xi,gi]=sc_selectg(Xi,genelist);    
%     si=s(c==i,:);
%     si=mean(si);
    %[Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
    %ctxt=Tct.C1_Cell_Type{1};
    
    ptsSelected=c==i;
    [Tct]=local_celltypebrushed(X,genelist,s,ptsSelected,...
          speciestag,organtag,databasetag);
    ctxt=Tct.C1_Cell_Type;

%     waitbar(1,f,'Finishing');
%     pause(1);
%     close(f);
    
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
            if size(s,2)>=3
                    %scatter3(s(ptsSelected,1),s(ptsSelected,2),s(ptsSelected,3),'x');
                    si=mean(s(ptsSelected,:));
                    text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
                         'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            elseif size(s,2)==2
                    %scatter(s(ptsSelected,1),s(ptsSelected,2),'x')                    
                    si=mean(s(ptsSelected,:));
                    text(si(:,1),si(:,2),sprintf('%s',ctxt),...
                         'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            end
            hold off
end
c_cell_type=string(cL(c));
export2wsdlg({'Save cell type list to variable named:'},...
    {'c_celltype'},{c_cell_type});
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
    [ax,bx]=view();
    [h]=i_gscatter3(s,c,methodid);
    title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
    view(ax,bx);
    [c,cL]=grp2idx(c);
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
    [Tct]=local_celltypebrushed(X,genelist,s,ptsSelected,...
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
            if size(s,2)>=3
                    scatter3(s(ptsSelected,1),s(ptsSelected,2),s(ptsSelected,3),'x');
                    si=mean(s(ptsSelected,:));
                    text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
                         'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            elseif size(s,2)==2
                    scatter(s(ptsSelected,1),s(ptsSelected,2),'x')                    
                    si=mean(s(ptsSelected,:));
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
    [markerlist]=sc_pickmarkers(X,genelist,1+ptsSelected,2);
    waitbar(1,f,'Finishing');
    pause(1);
    close(f);
    [ax,bx]=view();
    i_markergenespanel(X,genelist,s,...
        markerlist,numfig,9,ax,bx,...
        'Marker Genes for Selected Cells');
    pause(2);
    export2wsdlg({'Save marker list to variable named:'},...
            {'g_markerlist'},{markerlist});
end

function ShowMarkerGene(~,~)
    answer = questdlg('Select a gene to show expression?');
    if ~strcmp(answer,'Yes'), return; end
    gsorted=sort(genelist);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1        
        figure;
        sc_scattermarker(X,genelist,s,gsorted(indx),5);
%       [ax,bx]=view(); 
%       sc_markerscatter(X,genelist,gsorted(indx),s,3);
%       view(ax,bx);
    end
end

function ShowCellStats(~,~)
    answer = questdlg('Show cell states?');
    if ~strcmp(answer,'Yes'), return; end    
    [indx,tf] = listdlg('PromptString',{'Select statistics',...
    '',''},...    
    'SelectionMode','single',...
    'ListString',{'Library Size','Mt-reads Ratio',...
    'Mt-genes Expression','Cell Cycle Phase'});
    if tf==1        
        switch indx
            case 1
                ci=sum(X);
                ttxt="Library Size";
            case 2
                i=startsWith(genelist,'mt-','IgnoreCase',true);
                lbsz=sum(X,1);
                lbsz_mt=sum(X(i,:),1);
                ci=lbsz_mt./lbsz;
                ttxt="mtDNA%";
            case 3
                idx=startsWith(genelist,'mt-','IgnoreCase',true);
                n=sum(idx);
                if n>0
                    [ax,bx]=view();
                    if n<=9
                        i_markergenespanel(X,genelist,s,...
                            genelist(idx),[],9,ax,bx,'Mt-genes');
                    else
                        i_markergenespanel(X,genelist,s,...
                            genelist(idx),[],16,ax,bx,'Mt-genes');
                    end
                else
                    warndlg('No mt-genes found');
                end
                return;
            case 4   % "Cell Cycle Phase";
                if isempty(c_cell_cycle_phase)
                    f = waitbar(0,'Please wait...');
                    pause(.5)
                    waitbar(.67,f,'Processing your data');                
                    [cix]=run_cellcycle(X,genelist);
                    c_cell_cycle_phase=string(cix);
                    waitbar(1,f,'Finishing');
                    pause(1);
                    close(f);                    
                end                
                [ci,tx]=grp2idx(c_cell_cycle_phase);
                ttxt=sprintf('%s|',string(tx));
        end
            [ax,bx]=view();
            delete(h);
            h=i_gscatter3(s,ci);
            if indx==3, legend(tx); end
            view(ax,bx);
            title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)));
            hc=colorbar;
            hc.Label.String=ttxt;
            if indx==4
                labels = {'Save cell cycle phase to variable named:'};
                vars = {'c_cell_cycle_phase'};
                values = {c_cell_cycle_phase};
                export2wsdlg(labels,vars,values);
            end
    end
end

function SelectCellsByClass(~,~)
    answer = questdlg('Select cells by class?');
    if ~strcmp(answer,'Yes'), return; end    
    [indx,tf] = listdlg('PromptString',{'Select statistics','',''},...    
    'SelectionMode','single','ListString',...
     {'Custom input (C)','Cluster ID',...
            'Cell Type','Cell Cycle Phase'});
    if tf==1        
        [ax,bx]=view();
        switch indx
            case 1
                [indxx,tfx] = listdlg('PromptString',{'Select groups',...
                '',''},'SelectionMode','multiple','ListString',cL);
                if tfx==1
                    i=ismember(c,indxx);
                    [ax,bx]=view();                
                    sc_scatter(X(:,i),genelist,s(i,:),cL(c(i)));
                    view(ax,bx);
                end
            case 2
                if ~isempty(c_cluster_id)
                    [ci,cLi]=grp2idx(c_cluster_id);
                    [indxx,tfx] = listdlg('PromptString',{'Select groups',...
                    '',''},'SelectionMode','multiple','ListString',string(cLi));
                    if tfx==1
                        i=ismember(ci,indxx);
                        [ax,bx]=view();
                        sc_scatter(X(:,i),genelist,s(i,:),cLi(ci(i)));
                        view(ax,bx);
                    end
                else                    
                    errordlg('Class type undefined');
                    return;
                end
            case 3                
                if ~isempty(c_cell_type)
                    [ci,cLi]=grp2idx(c_cell_type);
                    [indxx,tfx] = listdlg('PromptString',{'Select groups',...
                    '',''},'SelectionMode','multiple','ListString',string(cLi));
                    if tfx==1                        
                        i=ismember(ci,indxx);                        
                        sc_scatter(X(:,i),genelist,s(i,:),cLi(ci(i)));
                        view(ax,bx);
                    end
                else
                    errordlg('Class type undefined');
                    return;
                end
            case 4
                if ~isempty(c_cell_cycle_phase)
                    [ci,cLi]=grp2idx(c_cell_cycle_phase);
                    [indxx,tfx] = listdlg('PromptString',{'Select groups',...
                    '',''},'SelectionMode','multiple','ListString',string(cLi));
                    if tfx==1                        
                        i=ismember(ci,indxx);                                       
                        sc_scatter(X(:,i),genelist,s(i,:),cLi(ci(i)));
                    end
                else
                    errordlg('Class type undefined');
                    return;
                end
        end
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
    X(:,ptsSelected)=[];
    s(ptsSelected,:)=[];
    c(ptsSelected)=[];    
    c_cell_idx(ptsSelected)=[];
    if ~isempty(c_cell_cycle_phase)
        c_cell_cycle_phase(ptsSelected)=[];
    end
    if ~isempty(c_cell_type)
        c_cell_type(ptsSelected)=[];
    end
    [ax,bx]=view();
    h=i_gscatter3(s,c);
    title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)));
    view(ax,bx);
end

function SaveX(~,~)
    answer = questdlg('Export & save data?');
    if ~strcmp(answer,'Yes'), return; end     
    labels = {'Save expression X to variable named:',...
              'Save embedding S to variable named:',...
              'Save group C to variable named:',...
              'Save cell index CELLIDX to variable named:'}; 
    vars = {'X_scatter','s_scatter','c_scatter','cellidx_scatter'};
    values = {X,s,c,c_cell_idx};
    export2wsdlg(labels,vars,values);
    %         assignin('base',sprintf('psexplorerT%d',...
    %                  psexplorer_timeid),t);
end

function DrawTrajectory(~,~)
        answer = questdlg('Which method?','Select Algorithm',...
            'splinefit (fast)','princurve (slow)',...
            'splinefit (fast)');
    if strcmp(answer,'splinefit (fast)')
        dim=1;
        [t,xyz1]=i_pseudotime_by_splinefit(s,dim,false);    
    else
        [t,xyz1]=i_pseudotime_by_princurve(s,false);
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
        values = {X, t, s};
        msgfig=export2wsdlg(labels,vars,values);
        %         assignin('base',sprintf('psexplorerT%d',...
        %                  psexplorer_timeid),t);
        uiwait(msgfig)
        answer = questdlg('View expression of selected genes', ...
            'Pseudotime Function', ...
            'Yes','No','Yes');
        switch answer
            case 'Yes'
                r=corr(t,X','type','spearman'); % Calculate linear correlation between gene expression profile and T
                [~,idxp]= maxk(r,4);  % Select top 4 positively correlated genes
                [~,idxn]= mink(r,3);  % Select top 3 negatively correlated genes
                selectedg=genelist([idxp idxn]);        
                figure;
                i_plot_pseudotimeseries(log2(X+1),...
                    genelist,t,selectedg);
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
    [t_mono,s_mono]=run_monocle(X);         waitbar(1,f,'Finishing');
    pause(1)
    close(f)   
    

        answer = questdlg('View Monocle DDRTree?', ...
            'Pseudotime View', ...
            'Yes','No','Yes');
        switch answer
            case 'Yes'
            [ax,bx]=view();
            cla(hAx);
            s=s_mono; c=t_mono;
            h=i_gscatter3(s,c);
            title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
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
    c_cluster_id=sc_cluster_s(s,k,'type',methodtag,'plotit',false);
    [c,cL]=grp2idx(c_cluster_id);
    waitbar(1,f,'Finishing');
    pause(1);
    close(f);
    [ax,bx]=view();
    h=i_gscatter3(s,c);
    view(ax,bx);
    title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
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
        [c_cluster_id]=sc_cluster_x(X,k,'type',methodtag);
    catch ME
        close(f);        
        errordlg(sprintf('%s: %s',...
            ME.identifier,ME.message));
        return;
    end
    [c,cL]=grp2idx(c_cluster_id);
    waitbar(1,f,'Finishing');
    pause(1);
    close(f);
    % hold off
    delete(h);
    h=i_gscatter3(s,c);
    title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
    
    labels = {'Save clusterid C to variable named:'}; 
    vars = {sprintf('c_clusterid_%s',methodtag)};
    values = {c_cluster_id};
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
            [ax,bx]=view();
            cla(hAx);
            h=i_gscatter3(s,c,methodid);
            title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
            view(ax,bx);
        else
            i_labelclusters;
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
    figure;
    for k=1:9
        if k<=max(c)
            subplot(3,3,k);
            i_gscatter3(s,c,3,cmv(idxx(k)));
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
                i_gscatter3(s,c,3,cmv(idxx(kk)));
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

function i_labelclusters
    stxtyes=false;
    if ~isempty(cL)
        answer = questdlg(sprintf('Label %d groups with index or text?',numel(cL)),...
            'Select Format','Index','Text','Text');
        if strcmp(answer,'Text')
            stxtyes=true; 
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
        si=s(c==i,:);
        si=mean(si);        
        if stxtyes
            stxt=sprintf('%s',cL{i});
        else
            stxt=sprintf('%d',i);
        end
        if size(s,2)==3
            text(si(:,1),si(:,2),si(:,3),stxt,...
                'fontsize',fsz,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
        else
            text(si(:,1),si(:,2),stxt,...
                'fontsize',fsz,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
        end
    end
    hold off
    % helpdlg(sprintf('%d clusters are labelled.',numel(cL)));
end

end


function [Tct]=local_celltypebrushed(X,genelist,s,...
                brushedData,species,organ,database)

if nargin<7, organ='panglaodb'; end
if nargin<6, organ='all'; end
if nargin<5, species='mouse'; end

if islogical(brushedData)
    i=brushedData;
else
    [~,i]=ismember(brushedData,s,'rows');
end
Xi=X(:,i);
[Xi,gi]=sc_selectg(Xi,genelist);
if strcmpi(database,'clustermole')
    %disp('Using clustermole marker database')
    [Tct]=sc_celltypecaller_new(Xi,gi,[],'species',species);
elseif strcmpi(database,'panglaodb')
    %disp('Using panglaodb marker database')
    [Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
end
end
