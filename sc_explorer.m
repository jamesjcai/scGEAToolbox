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
hFig = figure();
hAx = axes('Parent',hFig);
if size(s,2)>=5    
    i_view5d(s);
    hFig.Position(3)=hFig.Position(3)*2;
elseif size(s,2)==3
    scatter3(hAx,s(:,1),s(:,2),s(:,3),10);
elseif size(s,2)==2
    scatter(hAx,s(:,1),s(:,2),10);
end

tb = findall(hFig,'Type','uitoolbar'); % tb = uitoolbar(hFig);
pt = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_object01.gif'));
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Cell Type Explorer...';
pt.ClickedCallback = @MenuSelected1;

pt2 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_object01.gif'));
        map(9,:)=[0 0 0];
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Marker Gene Explorer...';
pt2.ClickedCallback = @MenuSelected2;

pt3 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_object02.gif'));
ptImage = ind2rgb(img,map);
pt3.CData = ptImage;
pt3.Tooltip = 'Pseudotime Explorer...';
pt3.ClickedCallback = @MenuSelected3;


       
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
 
%         cm = uicontextmenu(hFig);
%         m1 = uimenu(cm,'Text','Cell Type Explorer...');
%         m2 = uimenu(cm,'Text','Marker Gene Explorer...');
%         m3 = uimenu(cm,'Text','Pseudotime Explorer...');
%         % cm.ContextMenuOpeningFcn = @zzz;
%         hFig.ContextMenu = cm;        
%         m1.MenuSelectedFcn = @MenuSelected1;
%         m2.MenuSelectedFcn = @MenuSelected2;
%         m3.MenuSelectedFcn = @MenuSelected3;

add_3dcamera;

    function MenuSelected1(src,event)
        answer = questdlg('Which species?', ...
            'Select Species', ...
            'Mouse','Human','Mouse');
        switch answer
            case 'Human'
                speciesx="human";
            case 'Mouse'
                speciesx="mouse";
            otherwise
                return;
        end        
        answer = questdlg('Which algorithm?', ...
            'Select Method', ...
            'Alona','SingleR','Alona');
        switch answer
            case 'Alona'
                methodx="alona";
            case 'SingleR'
                methodx="singler";
            otherwise
                return;
        end        
        sc_celltypeexplorer(X,genelist,s,...
            'species',speciesx,"method",methodx);
    end

    function MenuSelected2(src,event)
        sc_markerexplorer(X,genelist,s);
    end
    function MenuSelected3(src,event)
        sc_pseudotimeexplorer(X,genelist,s);
    end
end



function i_view5d(s5d)
    if size(s5d,2)<5
        s5d=[s5d,zeros(size(s5d,1),5-size(s5d,2))];
    end
    i_view5dIDX=s5d(:,3);
    assignin('base','i_view5dIDX',i_view5dIDX);
    subplot(1,2,1);
    h1=scatter3(s5d(:,1),s5d(:,2),s5d(:,3),10);
    xlabel('dim 1')
    ylabel('dim 2')
    zlabel('dim 3')
    grid on
    box on
    h1.ZDataSource='i_view5dIDX';

    subplot(1,2,2);
    h2=scatter3(s5d(:,3),s5d(:,4),s5d(:,5),10);
    xlabel('dim 3')
    ylabel('dim 4')
    zlabel('dim 5')
    grid on
    box on
    h2.XDataSource='i_view5dIDX';
    hLD = linkdata('on');
    evalin('base','h=findobj(gcf,''type'',''axes'');');
    evalin('base','hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
    evalin('base','rotate3d on');
end
