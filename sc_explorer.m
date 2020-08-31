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
hFig = figure;
hAx = axes('Parent',hFig);

if size(s,2)>2
    scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
elseif size(s,2)==2
    scatter(hAx, s(:,1),s(:,2),10);
end

% defaultToolbar = findall(hFig,'Type','uitoolbar');
% pt = uipushtool(defaultToolbar);
% ptImage = rand(16,16,3);
% pt.CData = ptImage;
% pt.Tooltip = '___';
% pt.ClickedCallback = @showmenu;

    
        cm = uicontextmenu(hFig);
        m1 = uimenu(cm,'Text','Cell Type Explorer...');
        m2 = uimenu(cm,'Text','Marker Gene Explorer...');
        m3 = uimenu(cm,'Text','Pseudotime Explorer...');
        cm.ContextMenuOpeningFcn = @(src,event)disp('Context menu opened');
        hFig.ContextMenu = cm;        
        
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



