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
titxt=sprintf('Cell Type Explorer\n[species: %s; method: %s]',...
    species,method);
hFig = figure('Name','Cell Type Explorer');
hAx = axes('Parent',hFig);
if size(s,2)>=3
    scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
elseif size(s,2)==2
    scatter(hAx,s(:,1),s(:,2),10);
end
title(titxt);


tb = uitoolbar(hFig);
tt = uitoggletool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(mfilename('fullpath')),...
            '..','resources','tool_ellipse.gif'));
%  [img,map] = imread(fullfile(matlabroot,...
%              'toolbox','matlab','icons','tool_ellipse.gif'));
ptImage = ind2rgb(img,map);
tt.CData = ptImage;
tt.Tooltip = 'Click and then brush/select cells';
tt.ClickedCallback = @MenuSelected1;

    function MenuSelected1(src,~)
        state = src.State;        
        if strcmp(state,'on')
            hBr.Enable='on';
            % tt.CData = ptImage; % zeros(16,16,3);
        else
            hBr.Enable='off';
            % tt.CData = ptImage;
        end        
    end


hBr = brush(hFig);
% hBr.Enable='on';

hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,...
    species,organ,method};

gui.add_3dcamera(tb);

end


% ref: https://www.mathworks.com/matlabcentral/answers/385226-how-to-use-the-data-brush-tool-to-automatically-save-selected-points-in-multiple-line-plots

function onBrushAction(~,eventdata,X,genelist,s,species,organ,method)
%global ctexplorer_celltypeid
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
                    [Tct]=pkg.local_celltypebrushed(X,genelist,s,ptsSelected,species,organ);
                    ctxt=Tct.C1_Cell_Type;
                case 'singler'
                    % [~,i]=ismember(brushedData,s,'rows');
                    i=ptsSelected;
                    Xi=X(:,i);
                    [Xi,gi]=sc_selectg(Xi,genelist);
                    cx=run_singler(Xi,gi,species);
                    ctxt=unique(cx);
%                     Tct=tabulate(cx);
%                     %ii=grp2idx(ctxt);
%                     %mxii=mode(ii,'all');
%                     %ctxt=ctxt(unique(ii(mxii)));
%                     ctxt=unique(cx(mode(grp2idx(cx),'all')==grp2idx(cx)));
                end            
            [indx,tf] = listdlg('PromptString',{'Select cell type',...
            '',''},'SelectionMode','single','ListString',ctxt);
            if tf==1 
                ctxt=ctxt{indx};
            else
                return;
            end                
                
            %data = [hLines(k).XData(ptsSelected).' ...
            %    hLines(k).YData(ptsSelected).'];
            %assignin('base',names{k},data)
            
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
%            ctexplorer_celltypeid=ctexplorer_celltypeid+1;            
%             a=matlab.lang.makeValidName(ctxt);
%             a=extractBefore(a,min([10 strlength(a)]));
%             assignin('base',sprintf('ctexplorerX%d_%s',...
%                 ctexplorer_celltypeid,a),X(:,ptsSelected));
%             assignin('base',sprintf('ctexplorerT%d_%s',...
%                 ctexplorer_celltypeid,a),Tct);            
        end
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
% sc_celltypecaller
[Tct]=run.alona(Xi,gi,[],'species',species,'organ',organ);
end


