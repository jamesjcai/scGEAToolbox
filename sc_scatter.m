function sc_scatter(X,genelist,s,c,methodid)

if nargin<5, methodid=1; end
if nargin<4, c=ones(size(s,1),1); end

x=s(:,1);
y=s(:,2);

if iscell(c), c=grp2idx(c); end
kc=numel(unique(c));

if size(s,2)>=3, z=s(:,3); end

hFig = figure('Name','Pseudotime Explorer');
hAx = axes('Parent',hFig);


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
title(sprintf('%d x %d\n[genes x cells]',size(X,1),size(X,2)))
if kc<=5
    colormap(lines(kc));
else
    a=colormap('autumn');
    a(1,:)=[.8 .8 .8];
    colormap(a);
end
% add_3dcamera;

tb = uitoolbar(hFig);
pt = uipushtool(tb,'Separator','off');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-plot.gif'));
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Export & save data';
pt.ClickedCallback = @SaveX;


pt2 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-scatter.gif'));
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delete selected cells';
pt2.ClickedCallback = @DeleteSelectedCells;


pt3 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','list.gif'));         
ptImage = ind2rgb(img,map);
pt3.CData = ptImage;
pt3.Tooltip = 'Select a gene to show expression';
pt3.ClickedCallback = @showmkgene;

pt3a = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','list.gif'));         
ptImage = ind2rgb(img,map);
pt3a.CData = ptImage;
pt3a.Tooltip = 'Show cell states';
pt3a.ClickedCallback = @ShowCellstats;


pt4 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-stairs.gif'));
ptImage = ind2rgb(img,map);
pt4.CData = ptImage;
pt4.Tooltip = 'Marker genes of brushed cells';
pt4.ClickedCallback = @Brush4Markers;


pt5 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','brush.gif'));
ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Cell types of brushed cells';
pt5.ClickedCallback = @Brush4Celltypes;


add_3dcamera(tb);

% =========================
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
    f = waitbar(0,'Please wait...');
            pause(.5)
            waitbar(.67,f,'Processing your data');
            [markerlist]=sc_pickmarkers(X,genelist,1+ptsSelected,2);
            waitbar(1,f,'Finishing');
            pause(1)
            close(f)
            [ax,bx]=view();
             numfig=1;
             for kkk=1:numfig
            figure;
            for kk=1:min([9,length(markerlist)])
                subplot(3,3,kk)
                sc_markerscatter(X,genelist,...
                    markerlist(kk+9*(kkk-1)),s,3);
                view(ax,bx);   
            end
            end
%             mkexplorer_clustid=mkexplorer_clustid+1;
%             assignin('base',sprintf('mkexplorerL%d',...
%                 mkexplorer_clustid),markerlist);
end

function showmkgene(~,~)
    answer = questdlg('Select a gene to show expression?');
    if ~strcmp(answer,'Yes'), return; end
    gsorted=sort(genelist);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1        
        if size(s,1)==size(X,2)
        figure;
        sc_scattermarker(X,genelist,gsorted(indx),s,5);
%       [ax,bx]=view(); 
%       sc_markerscatter(X,genelist,gsorted(indx),s,3);
%       view(ax,bx);
        else
            errordlg('ERROR: size(s,1)!=size(X,2)')
        end
    end
end

function ShowCellstats(~,~)
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
              'Save group C to variable named:'}; 
    vars = {'X_scatter','c_scatter'};
    values = {X, c};
    msgfig=export2wsdlg(labels,vars,values);
    %         assignin('base',sprintf('psexplorerT%d',...
    %                  psexplorer_timeid),t);
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


