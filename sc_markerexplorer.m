function sc_markerexplorer(X,genelist,s,varargin)
%load cleandata.mat
%s=s_tsne;

   p = inputParser;
   addRequired(p,'X',@isnumeric);
   addRequired(p,'genelist',@isstring);
   addRequired(p,'s',@isnumeric);
   addOptional(p,'method',"ttest",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["ttest","mast"]));
   addOptional(p,'numfig',1,@isnumeric);
   parse(p,X,genelist,s,varargin{:});
   method=p.Results.method;
   numfig=p.Results.numfig;
   
  
global mkexplorer_clustid
mkexplorer_clustid=0;
hFig = figure('Name','Marker Gene Explorer');
hAx = axes('Parent',hFig);

if size(s,2)>=3
    scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
elseif size(s,2)==2
    scatter(hAx, s(:,1),s(:,2),10);
end

tb = uitoolbar(hFig);
tt = uitoggletool(tb);
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','tool_ellipse.gif'));
ptImage = ind2rgb(img,map);
tt.CData = ptImage;
tt.Tooltip = 'Click and then brush/select cells';
tt.ClickedCallback = @MenuSelected1;

    function MenuSelected1(src,event)
        state = src.State;        
        if strcmp(state,'on')
            hBr.Enable='on';
            tt.CData = zeros(16,16,3);
        else
            hBr.Enable='off';
            tt.CData = ptImage;
        end
    end


pt = uipushtool(tb,'Separator','off');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','profiler.gif'));
ptImage = ind2rgb(img,map);


% defaultToolbar = findall(hFig,'Type','uitoolbar');
% pt = uipushtool(defaultToolbar);
% ptImage = rand(16,16,3);
pt.CData = ptImage;
pt.Tooltip = 'Select a gene to show expression';
pt.ClickedCallback = @showmkgene;

function showmkgene(src,event)
    gsorted=sort(genelist);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1
        figure;
        sc_scattermarker(X,genelist,gsorted(indx),s,5);
%       [ax,bx]=view(); sc_markerscatter(X,genelist,gsorted(indx),s,3);
%       view(ax,bx);
    end
end

%h = uicontrol('Position',[5 5 150 30],'String','Calculate xdiff',...
%              'Callback', @JCal);

hBr = brush(hFig);
% hBr.Enable='on';
% hBr.Color = 'green';
hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,method,numfig};
end

% ref: https://www.mathworks.com/matlabcentral/answers/385226-how-to-use-the-data-brush-tool-to-automatically-save-selected-points-in-multiple-line-plots

function onBrushAction(~,eventdata,X,genelist,s,method,numfig)
global mkexplorer_clustid
% Extract plotted graphics objects
% Invert order because "Children" property is in reversed plotting order
hLines = flipud(eventdata.Axes.Children);

    % Loop through each graphics object
    for k = 1:1 %numel(hLines)
        % Check that the property is valid for that type of object
        % Also check if any points in that object are selected
        if isprop(hLines(k),'BrushData') && any(hLines(k).BrushData)
            % Output the selected data to the base workspace with assigned name
            ptsSelected = logical(hLines(k).BrushData.');
            switch method
                case 'mast'
                    [markerlist]=sc_pickmarkers2(X,genelist,1+ptsSelected,2);
                case 'ttest'
                    [markerlist]=sc_pickmarkers(X,genelist,1+ptsSelected,2);
            end
            [ax,bx]=view();
            for kkk=1:numfig
            figure;
            for kk=1:min([9,length(markerlist)])
                subplot(3,3,kk)
                sc_markerscatter(X,genelist,...
                    markerlist(kk+9*(kkk-1)),s,3);
                view(ax,bx);   
            end
            end
            mkexplorer_clustid=mkexplorer_clustid+1;
            assignin('base',sprintf('mkexplorerL%d',...
                mkexplorer_clustid),markerlist);
        end
    end
end



