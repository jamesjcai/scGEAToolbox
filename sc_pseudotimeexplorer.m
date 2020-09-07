function sc_pseudotimeexplorer(X,genelist,s,varargin)
%load cleandata.mat
%s=s_tsne;

   p = inputParser;
   addRequired(p,'X',@isnumeric);
   addRequired(p,'genelist',@isstring);
   addRequired(p,'s',@isnumeric);
   addOptional(p,'method',"splinefit",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["splinefit"]));
   addOptional(p,'dim',1,@isnumeric);
   parse(p,X,genelist,s,varargin{:});
   method=p.Results.method;
   dim=p.Results.dim;
   
  
global psexplorer_timeid
psexplorer_timeid=1;
hFig = figure('Name','Pseudotime Explorer');
hAx = axes('Parent',hFig);

if size(s,2)>=3
    hs=scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
elseif size(s,2)==2
    hs=scatter(hAx, s(:,1),s(:,2),10);
end

%defaultToolbar = findall(hFig,'Type','uitoolbar');
%pt = uipushtool(defaultToolbar);
tb = uitoolbar(hFig);
pt = uipushtool(tb,'Separator','off');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-plot.gif'));
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Plot pseudotime trajectory';
pt.ClickedCallback = @drawtrajectory;

pt2 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-scatter.gif'));
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delet selected cells';
pt2.ClickedCallback = @deleteselectedcells;
add_3dcamera(tb);

    function deleteselectedcells(~,~)
        data = hs.BrushData;
        ptsSelected=find(data);
        X(:,ptsSelected)=[];
        s(ptsSelected,:)=[];
        [a,b]=view();
        if size(s,2)>=3
            hs=scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
        elseif size(s,2)==2
            hs=scatter(hAx, s(:,1),s(:,2),10);
        end
        view(a,b);
    end

    function drawtrajectory(src,event)
        [t,xyz1]=i_pseudotime_by_splinefit(s,dim,false);
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
        assignin('base',sprintf('psexplorerT%d',...
                 psexplorer_timeid),t);
        r=corr(t,X','type','spearman'); % Calculate linear correlation between gene expression profile and T
        [~,idxp]= maxk(r,4);  % Select top 4 positively correlated genes
        [~,idxn]= mink(r,3);  % Select top 3 negatively correlated genes
        selectedg=genelist([idxp idxn]);
        % Plot expression profile of the 5 selected genes
        figure;
        i_plot_pseudotimeseries(log2(1+X),genelist,t,selectedg);
    end
end
