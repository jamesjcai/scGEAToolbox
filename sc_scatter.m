function sc_scatter(X,genelist,s,c,methodid)

if nargin<5, methodid=1; end
if nargin<4, c=ones(size(s,1),1); end
cellidx=(1:size(X,2))';

x=s(:,1);
y=s(:,2);
cL=[];
if iscell(c)||isstring(c), [c,cL]=grp2idx(c); end
kc=numel(unique(c));

if size(s,2)>=3, z=s(:,3); end

hFig = figure('Name','sc_scatter');
hAx = axes('Parent',hFig);
[h]=i_gscatter3(s,c,methodid);
%{
switch methodid
    case 1
        if size(s,2)==2
           h=scatter(hAx,x,y,10,c);
        elseif size(s,2)>=3
           h=scatter3(hAx,x,y,z,10,c);
        end
    case 2
        if size(s,2)==2
            h=gscatter(x,y,c,[],[],10);
        elseif size(s,2)>=3
            h=gscatter3b(x,y,z,c);
        end
end
%}
title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
if kc<=20
    colormap(lines(kc));
else
    a=colormap('autumn');
    a(1,:)=[.8 .8 .8];
    colormap(a);
end
% add_3dcamera;

dt=datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};


tb = uitoolbar(hFig);
pt3 = uipushtool(tb,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','list.gif'));         
ptImage = ind2rgb(img,map);
pt3.CData = ptImage;
pt3.Tooltip = 'Select a gene to show expression';
pt3.ClickedCallback = @ShowMarkerGene;

pt3a = uipushtool(tb,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','list2.gif'));         
ptImage = ind2rgb(img,map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Show cell states';
pt3a.ClickedCallback = @ShowCellStats;

% ------------------

ptlabelclusters = uitoggletool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
             'toolbox','matlab','icons','plotpicker-scatter.gif'));
% map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
ptImage = ind2rgb(img,map);
ptlabelclusters.CData = ptImage;
ptlabelclusters.Tooltip = 'Label clusters';
ptlabelclusters.ClickedCallback = @LabelClusters;


ptaddcluster = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','plotpicker-pzmap.gif'));
ptImage = ind2rgb(img,map);
ptaddcluster.CData = ptImage;
ptaddcluster.Tooltip = 'Add brushed cells to a new cluster';
ptaddcluster.ClickedCallback = @Brushed2Cluster;


ptcluster = uipushtool(tb,'Separator','off');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','plotpicker-dendrogram.gif'));
ptImage = ind2rgb(img,map);
ptcluster.CData = ptImage;
ptcluster.Tooltip = 'Clustering cells';
ptcluster.ClickedCallback = @ClusterCells;


ptclustertype = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
             'toolbox','matlab','icons','plotpicker-contour.gif'));
ptImage = ind2rgb(img,map);
ptclustertype.CData = ptImage;
ptclustertype.Tooltip = 'Cell types of clusters';
ptclustertype.ClickedCallback = @CellTypeClusters;

pt5 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','brush.gif'));
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Cell types of brushed cells';
pt5.ClickedCallback = @Brush4Celltypes;


pt4 = uipushtool(tb,'Separator','on');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-stairs.gif'));
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','plotpicker-kagi.gif'));
ptImage = ind2rgb(img,map);
pt4.CData = ptImage;
pt4.Tooltip = 'Marker genes of brushed cells';
pt4.ClickedCallback = @Brush4Markers;


ptpseudotime = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','plotpicker-comet.gif'));
ptImage = ind2rgb(img,map);
ptpseudotime.CData = ptImage;
ptpseudotime.Tooltip = 'Plot pseudotime trajectory';
ptpseudotime.ClickedCallback = @drawtrajectory;


pt2 = uipushtool(tb,'Separator','on');
%[img,map] = imread(fullfile(matlabroot,...
%            'toolbox','matlab','icons','plotpicker-scatter.gif'));
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','plotpicker-qqplot.gif'));  
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delete selected cells';
pt2.ClickedCallback = @DeleteSelectedCells;

pt = uipushtool(tb,'Separator','off');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','plotpicker-plot.gif'));
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','export.gif'));         
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Export & save data';
pt.ClickedCallback = @SaveX;


pt5 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','plotpicker-geobubble.gif'));
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Refresh';
pt5.ClickedCallback = @RefreshAll;


add_3dcamera(tb);

%warning off
%WinOnTop(hFig,true);
%warning on

% =========================
function RefreshAll(~,~)
    hold off
    h=i_gscatter3(s,c,methodid);
    title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
    ptlabelclusters.State='off';
    tb.Visible='off';
    tb.Visible='on';
end

function CellTypeClusters(~,~)
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
    
    prompt = {'Enter number of panels {1,2..10}'};
    dlgtitle = 'Panel of 9 Genes';
    answer = inputdlg(prompt,dlgtitle,[1 40],{'1'});
    if isempty(answer)
        return;
    else
        numfig=str2double(answer{1});
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
            pause(1)
            close(f)
            [ax,bx]=view();
            % numfig=1;
            for kkk=1:numfig
                figure;
                for kk=1:min([9,length(markerlist)])
                    subplot(3,3,kk)                    
                    sc_scattermarker(X,genelist,s,...
                        markerlist(kk+9*(kkk-1)),3,5,false);
                    view(ax,bx);   
                end
            end
    export2wsdlg({'Save marker list to variable named:'},...
        {'g_markerlist'},{markerlist});
    
%             mkexplorer_clustid=mkexplorer_clustid+1;
%             assignin('base',sprintf('mkexplorerL%d',...
%                 mkexplorer_clustid),markerlist);
end

function ShowMarkerGene(~,~)
    answer = questdlg('Select a gene to show expression?');
    if ~strcmp(answer,'Yes'), return; end
    gsorted=sort(genelist);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1        
        if size(s,1)==size(X,2)
        figure;
        sc_scattermarker(X,genelist,s,gsorted(indx),5);
%       [ax,bx]=view(); 
%       sc_markerscatter(X,genelist,gsorted(indx),s,3);
%       view(ax,bx);
        else
            errordlg('ERROR: size(s,1)!=size(X,2)')
        end
    end
end

function ShowCellStats(~,~)
    answer = questdlg('Show cell states?');
    if ~strcmp(answer,'Yes'), return; end    
    [indx,tf] = listdlg('PromptString',{'Select statistics',...
    '',''},...    
    'SelectionMode','single',...
    'ListString',{'Library Size','Mt-reads Ratio','Cell Cycle Phase'});
    if tf==1
        if size(s,1)==size(X,2)
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
                % ttxt="Cell Cycle Phase";
                f = waitbar(0,'Please wait...');
                pause(.5)
                waitbar(.67,f,'Processing your data');                
                [cix]=run_cellcycle(X,genelist);
                waitbar(1,f,'Finishing');
                pause(1)
                close(f)
                [ci,tx]=grp2idx(cix);
                ttxt=sprintf('%s|',string(tx));
        end
            
            [ax,bx]=view();
                figure;
                if size(s,2)==2
                    scatter(s(:,1),s(:,2),5,ci,'filled');
                else
                    scatter3(s(:,1),s(:,2),s(:,3),5,ci,'filled');
                end
                title(ttxt);
%               axx=colormap('autumn');
%               % axx(1,:)=[.8 .8 .8];
%               colormap(axx);
                hc=colorbar;
                hc.Label.String=ttxt;
                % ylabel(hc,ttxt,'Rotation',270)
                view(ax,bx);
                
                if indx==3
                    colormap(lines(3)); 
                    pause(2);
                    labels = {'Save cell cycle phase to variable named:'}; 
                    vars = {'c_cell_cycle_phase'};
                    values = {cix};
                    msgfig=export2wsdlg(labels,vars,values);
                end                
        else
            errordlg('ERROR: size(s,1)!=size(X,2)')
        end
    end
end

function DeleteSelectedCells(~,~)
    answer = questdlg('Delete selected cells?');
    if ~strcmp(answer,'Yes'), return; end 
    ptsSelected = logical(h.BrushData.');    
    if ~any(ptsSelected)
        warndlg("No cells are selected.");
        return;
    end    
    data = h.BrushData;
    ptsSelected=find(data);
    X(:,ptsSelected)=[];
    s(ptsSelected,:)=[];
    c(ptsSelected)=[];
    cellidx(ptsSelected)=[];
    
    [a,b]=view();
    if size(s,2)>=3
        h=scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
    elseif size(s,2)==2
        h=scatter(hAx, s(:,1),s(:,2),10);
    end
    view(a,b);
end

function SaveX(~,~)
    answer = questdlg('Export & save data?');
    if ~strcmp(answer,'Yes'), return; end     
    labels = {'Save expression X to variable named:',...
              'Save embedding S to variable named:',...
              'Save group C to variable named:',...
              'Save cell index CELLIDX to variable named:'}; 
    vars = {'X_scatter','s_scatter','c_scatter','cellidx_scatter'};
    values = {X,s,c,cellidx};
    msgfig=export2wsdlg(labels,vars,values);
    %         assignin('base',sprintf('psexplorerT%d',...
    %                  psexplorer_timeid),t);
end

function drawtrajectory(~,~)
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

function ClusterCells(~,~)
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
    hold off
    c=sc_clustshow(s,k,'type',methodtag,'plotit',false);
    h=i_gscatter3(s,c);
    title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
    answer = questdlg('Label clusters?');
    if strcmp(answer,'Yes')
        i_labelclusters;
    end    
end

function LabelClusters(src,~)
        state = src.State;
        if strcmp(state,'off')
            [ax,bx]=view();
            hold off
            h=i_gscatter3(s,c,methodid);
            title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
            view(ax,bx);
        else
            i_labelclusters;
        end
end

function [txt]=i_myupdatefcnx(~,event_obj)
    % pos = event_obj.Position;
    idx = event_obj.DataIndex;
    txt =sprintf('class=%d',c(idx));
end

function i_labelclusters
    stxtyes=false;
    if ~isempty(cL)
        answer = questdlg('Label with index or text?','Select Format','Index','Text','Text');
        if strcmp(answer,'Text'), stxtyes=true; end
    end    
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
                'fontsize',20,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
        else
            text(si(:,1),si(:,2),stxt,...
                'fontsize',20,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
        end
    end
    hold off
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
    disp('Using clustermole marker database')
    [Tct]=sc_celltypecaller_new(Xi,gi,[],'species',species);
elseif strcmpi(database,'panglaodb')
    disp('Using panglaodb marker database')
    [Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
end
end



