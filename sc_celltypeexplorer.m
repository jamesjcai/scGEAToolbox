function sc_celltypeexplorer(X,genelist,s,varargin)
%load cleandata.mat
%s=s_tsne;

   p = inputParser;
   addRequired(p,'X',@isnumeric);
   addRequired(p,'genelist',@isstring);
   addRequired(p,'s',@isnumeric);
   addOptional(p,'species',"mouse",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["human","mouse"]));
   addOptional(p,'organ',"all",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["all","heart","immunesystem","brain","pancreas"]));
   addOptional(p,'method',"alona",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["alona","singler"]));
   parse(p,X,genelist,s,varargin{:});
   species=p.Results.species;
   organ=p.Results.organ;
   method=p.Results.method;
   
   if strcmpi(species,'mm')
       species="mouse";
   elseif strcmpi(species,'hs')
       species="human";
   end
   
global ctexplorer_celltypeid
ctexplorer_celltypeid=0;
hFig = figure;
hAx = axes('Parent',hFig);
scatter3(hAx, s(:,1),s(:,2),s(:,3),10);

hBr = brush(hFig);
% hBr.Enable='on';
hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,species,organ,method};
end


% ref: https://www.mathworks.com/matlabcentral/answers/385226-how-to-use-the-data-brush-tool-to-automatically-save-selected-points-in-multiple-line-plots

function onBrushAction(~,eventdata,X,genelist,s,species,organ,method)
global ctexplorer_celltypeid
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
            % find(ptsSelected)
            
                switch lower(method)
                case 'alona'
                    %[Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
                    [Tct]=sc_celltypebrushed(X,genelist,s,ptsSelected,species,organ);
                    ctxt=Tct.C1_Cell_Type{1};
                case 'singler'
                    % [~,i]=ismember(brushedData,s,'rows');
                    i=ptsSelected;
                    Xi=X(:,i);
                    [Xi,gi]=sc_selectg(Xi,genelist);
                    cx=run_singler(Xi,gi,species);
                    Tct=tabulate(cx);
                    %ii=grp2idx(ctxt);
                    %mxii=mode(ii,'all');
                    %ctxt=ctxt(unique(ii(mxii)));
                    ctxt=unique(cx(mode(grp2idx(cx),'all')==grp2idx(cx)));
                end    
            
            %data = [hLines(k).XData(ptsSelected).' ...
            %    hLines(k).YData(ptsSelected).'];
            %assignin('base',names{k},data)
            hold on
            scatter3(s(ptsSelected,1),s(ptsSelected,2),s(ptsSelected,3),'x')
            si=mean(s(ptsSelected,:));
            text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
                 'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
            hold off
            
            a=matlab.lang.makeValidName(ctxt);
            a=extractBefore(a,min([10 strlength(a)]));
            ctexplorer_celltypeid=ctexplorer_celltypeid+1;
            assignin('base',sprintf('ctexplorerX%d_%s',...
                ctexplorer_celltypeid,a),X(:,ptsSelected));
            assignin('base',sprintf('ctexplorerT%d_%s',...
                ctexplorer_celltypeid,a),Tct);            
        end
    end
end



