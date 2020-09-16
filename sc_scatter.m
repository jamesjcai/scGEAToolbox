function sc_scatter(X,genelist,s,c,methodid)

if nargin<5, methodid=1; end
if nargin<4, c=ones(size(s,1),1); end

x=s(:,1);
y=s(:,2);

if iscell(c), c=grp2idx(c); end
kc=numel(unique(c));

if size(s,2)>=3, z=s(:,3); end

hFig = figure('Name','Pseudotime Explorer');
hAx = axes('Parent',hFig);


switch methodid
    case 1
        if size(s,2)==2
           h=scatter(hAx,x,y,10,c);
        elseif size(s,2)>=3
           h=scatter3(hAx,x,y,z,10,c);
        end
    case 2
        if size(s,2)==2
            h=gscatter(x,y,c,[],[],10);
        elseif size(s,2)>=3
            h=gscatter3b(x,y,z,c);
        end
end

if kc<=5
    colormap(lines(kc));
else
    a=colormap('autumn');
    a(1,:)=[.8 .8 .8];
    colormap(a);
end
% add_3dcamera;

tb = uitoolbar(hFig);
pt = uipushtool(tb,'Separator','off');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-plot.gif'));
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Save X';
pt.ClickedCallback = @saveX;

pt2 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','plotpicker-scatter.gif'));
ptImage = ind2rgb(img,map);
pt2.CData = ptImage;
pt2.Tooltip = 'Delet selected cells';
pt2.ClickedCallback = @deleteselectedcells;


pt3 = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'private','list.gif'));         
ptImage = ind2rgb(img,map);
pt3.CData = ptImage;
pt3.Tooltip = 'Select a gene to show expression';
pt3.ClickedCallback = @showmkgene;

add_3dcamera(tb);




function showmkgene(~,~)
    gsorted=sort(genelist);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1        
        if size(s,1)==size(X,2)
        figure;
        sc_scattermarker(X,genelist,gsorted(indx),s,5);
%       [ax,bx]=view(); 
%       sc_markerscatter(X,genelist,gsorted(indx),s,3);
%       view(ax,bx);
        else
            errordlg('ERROR: size(s,1)!=size(X,2)')
        end
    end
end




    function deleteselectedcells(~,~)
        data = h.BrushData;
        ptsSelected=find(data);
        X(:,ptsSelected)=[];
        s(ptsSelected,:)=[];
        c(ptsSelected)=[];
        [a,b]=view();
        if size(s,2)>=3
            h=scatter3(hAx, s(:,1),s(:,2),s(:,3),10);
        elseif size(s,2)==2
            h=scatter(hAx, s(:,1),s(:,2),10);
        end
        view(a,b);
    end

    function saveX(~,~)
        labels = {'Save expression X to variable named:',...
                  'Save group C to variable named:'}; 
        vars = {'X_scatter','c_scatter'};
        values = {X, c};
        msgfig=export2wsdlg(labels,vars,values);
        %         assignin('base',sprintf('psexplorerT%d',...
        %                  psexplorer_timeid),t);
       
    end

end


