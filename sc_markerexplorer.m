function sc_markerexplorer(X,genelist,s,varargin)
%load cleandata.mat
%s=s_tsne;

   p = inputParser;
   addRequired(p,'X',@isnumeric);
   addRequired(p,'genelist',@isstring);
   addRequired(p,'s',@isnumeric);
   addOptional(p,'method',"ttest",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["ttest","mast"]));
   parse(p,X,genelist,s,varargin{:});
   method=p.Results.method;
   
  
global mkexplorer_clustid
mkexplorer_clustid=0;
hFig = figure;
hAx = axes('Parent',hFig);
scatter3(hAx, s(:,1),s(:,2),s(:,3),10);

hBr = brush(hFig);
% hBr.Enable='on';
hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,method};
end


% ref: https://www.mathworks.com/matlabcentral/answers/385226-how-to-use-the-data-brush-tool-to-automatically-save-selected-points-in-multiple-line-plots

function onBrushAction(~,eventdata,X,genelist,s,method)
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
            figure;
            for kk=1:min([9,length(markerlist)])
                subplot(3,3,kk)
                sc_scattermarker(X,genelist,markerlist(kk),s,3);
                view(ax,bx);            
            end
            assignin('base',sprintf('mkexplorerL%d',...
                mkexplorer_clustid),markerlist);
        end
    end
end



