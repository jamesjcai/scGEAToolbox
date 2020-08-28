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

if size(s,2)==3
    scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
elseif size(s,2)==2
    scatter(hAx, s(:,1),s(:,2),10);
end

defaultToolbar = findall(hFig,'Type','uitoolbar');
pt = uipushtool(defaultToolbar);
ptImage = rand(16,16,3);
pt.CData = ptImage;
pt.Tooltip = '___';
% pt.ClickedCallback = @showmkgene;

%hBr = brush(hFig);
%hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,method,numfig};

end

