function sc_explorer(X,genelist,s,varargin)
%load cleandata.mat
%s=s_tsne;

   p = inputParser;
   addRequired(p,'X',@isnumeric);
   addRequired(p,'genelist',@isstring);
   addRequired(p,'s',@isnumeric);
   addOptional(p,'method',"marker",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["marker","celltype","pseudotime"]));
   addOptional(p,'numfig',1,@isnumeric);
   parse(p,X,genelist,s,varargin{:});
   method=p.Results.method;
   numfig=p.Results.numfig;   
  
global mkexplorer_clustid
mkexplorer_clustid=0;

if size(s,2)>=5
    [hFig]=sc_explorer5d(s,X,genelist);    
elseif size(s,2)==3
    hFig = figure;
    hAx = axes('Parent',hFig);
    scatter3(hAx,s(:,1),s(:,2),s(:,3),10);
elseif size(s,2)==2
    hFig = figure;
    hAx = axes('Parent',hFig);
    scatter(hAx, s(:,1),s(:,2),10);
end

%defaultToolbar = findall(hFig,'Type','uitoolbar');
tb = uitoolbar(hFig);
pt = uipushtool(tb);
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_object01.gif'));
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Cell Type Explorer...';
pt.ClickedCallback = @MenuSelected1;

pt2 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_object02.gif'));
        map(9,:)=[0 0 0];
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Marker Gene Explorer...';
pt2.ClickedCallback = @MenuSelected2;

pt3 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_grid.gif'));
ptImage = ind2rgb(img,map);
pt3.CData = ptImage;
pt3.Tooltip = 'Pseudotime Explorer...';
pt3.ClickedCallback = @MenuSelected3;



        cm = uicontextmenu(hFig);
        m1 = uimenu(cm,'Text','Cell Type Explorer...');
        m2 = uimenu(cm,'Text','Marker Gene Explorer...');
        m3 = uimenu(cm,'Text','Pseudotime Explorer...');
        % cm.ContextMenuOpeningFcn = @zzz;
        hFig.ContextMenu = cm;        
m1.MenuSelectedFcn = @MenuSelected1;
m2.MenuSelectedFcn = @MenuSelected2;
m3.MenuSelectedFcn = @MenuSelected3;

       
%hBr = brush(hFig);
%hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,method,numfig};

m = uimenu('Text','&Explorers');
mitem1 = uimenu(m,'Text','&Cell Type Explorer...');
mitem2 = uimenu(m,'Text','&Marker Gene Explorer...');
mitem3 = uimenu(m,'Text','&Pseudotime Explorer...');
mitem1.Accelerator = 'C';
mitem2.Accelerator = 'M';
mitem3.Accelerator = 'P';
mitem1.MenuSelectedFcn = @MenuSelected1;
mitem2.MenuSelectedFcn = @MenuSelected2;
mitem3.MenuSelectedFcn = @MenuSelected3;
 
    function MenuSelected1(src,event)
        sc_celltypeexplorer(X,genelist,s);
    end
    function MenuSelected2(src,event)
        sc_markerexplorer(X,genelist,s);
    end
    function MenuSelected3(src,event)
        sc_pseudotimeexplorer(X,genelist,s);
    end
end



