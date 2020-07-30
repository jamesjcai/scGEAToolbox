function sc_celltypeexplorer(X,genelist,s,varargin)
%load cleandata.mat
%s=s_tsne;

   p = inputParser;
   addOptional(p,'species',"human",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["human","mouse"]));
   parse(p,varargin{:});
   species=p.Results.species;
   
global ctexplorer_celltypeid
ctexplorer_celltypeid=0;
hFig = figure;
hAx = axes('Parent',hFig);
scatter3(hAx, s(:,1),s(:,2),s(:,3),10);

hBr = brush(hFig);
% hBr.Enable='on';
hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,species};
end


% ref: https://www.mathworks.com/matlabcentral/answers/385226-how-to-use-the-data-brush-tool-to-automatically-save-selected-points-in-multiple-line-plots

function onBrushAction(~,eventdata,X,genelist,s,species)
global ctexplorer_celltypeid
% Extract plotted graphics objects
% Invert order because "Children" property is in reversed plotting order
hLines = flipud(eventdata.Axes.Children);

    % Loop through each graphics object
    for k = 1:numel(hLines)
        % Check that the property is valid for that type of object
        % Also check if any points in that object are selected
        if isprop(hLines(k),'BrushData') && any(hLines(k).BrushData)
            % Output the selected data to the base workspace with assigned name
            ptsSelected = logical(hLines(k).BrushData.');
            % find(ptsSelected)
            [Tct]=sc_celltypebrushed(X,genelist,s,ptsSelected,species);
            %data = [hLines(k).XData(ptsSelected).' ...
            %    hLines(k).YData(ptsSelected).'];
            %assignin('base',names{k},data)
            hold on
            scatter3(s(ptsSelected,1),s(ptsSelected,2),s(ptsSelected,3),'x')
            si=mean(s(ptsSelected,:));
            text(si(:,1),si(:,2),si(:,3),sprintf('%s',Tct.C1_Cell_Type{1}),...
                 'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            hold off
            
            a=matlab.lang.makeValidName(Tct.C1_Cell_Type{1});
            a=extractBefore(a,min([10 length(a)]));
            ctexplorer_celltypeid=ctexplorer_celltypeid+1;
            assignin('base',sprintf('ctexplorerX%d_%s',...
                ctexplorer_celltypeid,a),X(:,ptsSelected));
            assignin('base',sprintf('ctexplorerT%d_%s',...
                ctexplorer_celltypeid,a),Tct);            
        end
    end
end

