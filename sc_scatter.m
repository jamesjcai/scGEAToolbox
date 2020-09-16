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
pt.Tooltip = 'Save X';
pt.ClickedCallback = @saveX;


pt2 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-scatter.gif'));
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delet selected cells';
pt2.ClickedCallback = @deleteselectedcells;


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
pt3a.Tooltip = 'Show cell stats';
pt3a.ClickedCallback = @ShowCellstats;


pt4 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-stairs.gif'));
ptImage = ind2rgb(img,map);
pt4.CData = ptImage;
pt4.Tooltip = 'Marker gene brush';
pt4.ClickedCallback = @Brush4Markers;


pt5 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','brush.gif'));ptImage = ind2rgb(img,map);
pt5.CData = ptImage;
pt5.Tooltip = 'Cell type brush';
pt5.ClickedCallback = @Brush4Celltypes;


add_3dcamera(tb);

% =========================

function Brush4Celltypes(~,~)
            f = waitbar(0,'Please wait...');
            pause(.5)
            waitbar(.67,f,'Processing your data');
            ptsSelected = logical(h.BrushData.');
            if ~any(ptsSelected)
                waitbar(1,f,'Finishing');
                close(f)
                warndlg("No cells are selected.");
                return;
            end
            species="mouse";
            organ="all";
            [Tct]=local_celltypebrushed(X,genelist,s,ptsSelected,species,organ);
            ctxt=Tct.C1_Cell_Type{1};            
            [markerlist]=sc_pickmarkers(X,genelist,1+ptsSelected,2);
            waitbar(1,f,'Finishing');
            pause(1)
            close(f)
            
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
            f = waitbar(0,'Please wait...');
            pause(.5)
            waitbar(.67,f,'Processing your data');
            ptsSelected = logical(h.BrushData.');
            if ~any(ptsSelected)
                waitbar(1,f,'Finishing');
                close(f)
                warndlg("No cells are selected.");
                return;
            end            
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
    [indx,tf] = listdlg('PromptString',{'Select statistics',...
    '',''},...    
    'SelectionMode','single',...
    'ListString',{'Library Size','Mt-reads Ratio'});
    if tf==1
        if size(s,1)==size(X,2)
        switch indx
            case 1
                ci=sum(X);
            case 2
                i=startsWith(genelist,'mt-','IgnoreCase',true);
                lbsz=sum(X,1);
                lbsz_mt=sum(X(i,:),1);
                ci=lbsz_mt./lbsz;
        end
            
            [ax,bx]=view();
                figure;
                if size(s,2)==2
                    scatter(s(:,1),s(:,2),[],ci,'filled');
                else
                    scatter3(s(:,1),s(:,2),s(:,3),[],ci,'filled');
                end        
                axx=colormap('autumn');
                % axx(1,:)=[.8 .8 .8];
                colormap(axx);
                colorbar;
            view(ax,bx);
        else
            errordlg('ERROR: size(s,1)!=size(X,2)')
        end
    end
end


function deleteselectedcells(~,~)
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

function saveX(~,~)
    labels = {'Save expression X to variable named:',...
              'Save group C to variable named:'}; 
    vars = {'X_scatter','c_scatter'};
    values = {X, c};
    msgfig=export2wsdlg(labels,vars,values);
    %         assignin('base',sprintf('psexplorerT%d',...
    %                  psexplorer_timeid),t);
end

end



function [Tct]=local_celltypebrushed(X,genelist,s,brushedData,species,organ)

% USAGE:

% s=sc_tsne(X,3);
% figure; sc_cellscatter(s)
% % get brushedData
% [Tct]=sc_celltypesbrushed(X,genelist,s,brushedData)
if nargin<6, organ='all'; end
if nargin<5, species='human'; end

if islogical(brushedData)
    i=brushedData;
else
    [~,i]=ismember(brushedData,s,'rows');
end
Xi=X(:,i);
[Xi,gi]=sc_selectg(Xi,genelist);
[Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
end


