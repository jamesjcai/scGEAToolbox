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

% ===========================================
%{
bg = uibuttongroup(hFig, 'Visible','off',...
                  'Units','pixels',...
                  'Position',[5 5 70 72],...
                  'SelectionChangedFcn',@bselection,...
                  'title','Species');
              
% Create three radio buttons in the button group.
r1 = uicontrol(bg,'Style',...
                  'radiobutton',...
                  'String','Mouse',...
                  'Position',[5 25 70 25],...
                  'HandleVisibility','off');
              
r2 = uicontrol(bg,'Style','radiobutton',...
                  'String','Human',...
                  'Position',[5 5 70 25],...
                  'HandleVisibility','off');
bg.Visible = 'on';

    function bselection(source,event)
       % disp(['Previous: ' event.OldValue.String]);
       disp(['Current: ' event.NewValue.String]);
       disp('------------------');
       if strcmp(event.NewValue.String,"Mouse")
           species="mouse";
       else
           species="human";
       end
    end

%}
% =============

%{
bg2 = uibuttongroup(hFig, 'Visible','off',...
                  'Units','pixels',...
                  'Position',[5 75 70 72],...
                  'SelectionChangedFcn',@bselection2,...
                  'title','Method');
              
r1x = uicontrol(bg2,'Style',...
                  'radiobutton',...
                  'String','Alona',...
                  'Position',[5 25 70 25],...
                  'HandleVisibility','off');
              
r2x = uicontrol(bg2,'Style','radiobutton',...
                  'String','SingleR',...
                  'Position',[5 5 70 25],...
                  'HandleVisibility','off');
bg2.Visible = 'on';

    function bselection2(~,event2)
       % disp(['Previous: ' event.OldValue.String]);
       disp(['Current: ' event2.NewValue.String]);
       disp('------------------');
       if strcmp(event2.NewValue.String,"Alona")
           method="alona";
       else
           method="singler";
       end
    end
%}

% ===========================================

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


hBr = brush(hFig);
% hBr.Enable='on';

hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,...
    species,organ,method};

add_3dcamera(tb);

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
                    %[Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);                    
                    [Tct]=local_celltypebrushed(X,genelist,s,ptsSelected,species,organ);
                    ctxt=Tct.C1_Cell_Type{1};
                case 'singler'
                    % [~,i]=ismember(brushedData,s,'rows');
                    disp('xxxx');
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
[Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
end


