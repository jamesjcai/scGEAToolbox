function hLink=sc_explorer3(s1,X,genelist)
if nargin<3, genelist=[]; end
if nargin<2, X=[]; end
if size(s1,2)<5
    s1=[s1,zeros(size(s1,1),5-size(s1,2))];
end

% sc_explorer3(s_tsne_5d)

explorer2IDX=s1(:,3);
assignin('base','explorer2IDX',explorer2IDX);


hFig=figure;
%https://www.mathworks.com/matlabcentral/answers/153-if-i-have-two-plots-on-the-same-figure-window-how-do-i-use-the-brush-tool-to-highlight-one-data-poi
%https://www.mathworks.com/matlabcentral/answers/385300-how-to-set-the-datasource-of-a-histogram-programmatically

ax1 = subplot(1,2,1);
h1=scatter3(s1(:,1),s1(:,2),s1(:,3),10);
xlabel('dim 1')
ylabel('dim 2')
zlabel('dim 3')
grid on
box on
h1.ZDataSource='explorer2IDX';

ax2 = subplot(1,2,2);
h2=scatter3(s1(:,3),s1(:,4),s1(:,5),10);
xlabel('dim 3')
ylabel('dim 4')
zlabel('dim 5')
grid on
box on
h2.XDataSource='explorer2IDX';
hLD = linkdata('on');


% figure
% ax1 = subplot(2,1,1);
% [X1,Y1,Z1] = peaks;
% surf(X1,Y1,Z1)
% 
% ax2 = subplot(2,1,2);
% [X2,Y2,Z2] = peaks(10);
% surf(X2,Y2,Z2)

%assignin('base','ax1',ax1);
%assignin('base','ax2',ax2);


%hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); 
%rotate3d on

%h=findobj(gcf,'type','axes');
%hlink = linkprop(h,{'CameraPosition','CameraUpVector'}); 
%rotate3d on

evalin('base','h=findobj(gcf,''type'',''axes'');');
evalin('base','hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
evalin('base','rotate3d on');
hFig.Position(3)=hFig.Position(3)*2;


if ~isempty(X)&& ~isempty(genelist)
tb = uitoolbar(hFig);
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
end

function showmkgene(src,event)
    gsorted=sort(genelist);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1  
        [ax,bx]=view();
        figure;
        sc_markerscatter(X,genelist,gsorted(indx),s1,3);
        view(ax,bx);  
    end
end


tt = uitoggletool(tb,'State','on','Separator','on');
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','greencircleicon.gif'));
ptImage = ind2rgb(img,map);
tt.CData = ptImage;
tt.Tooltip = 'Link panels';
tt.ClickedCallback = @MenuSelected1;

    function MenuSelected1(src,event)
        state = src.State;        
        if strcmp(state,'on')            
            evalin('base','hlink.Enabled=''on''');
            [img,map] = imread(fullfile(matlabroot,...
                        'toolbox','matlab','icons','greencircleicon.gif'));
            ptImage = ind2rgb(img,map);            
            tt.CData = ptImage;
        else            
            [img,map] = imread(fullfile(matlabroot,...
                        'toolbox','matlab','icons','tool_ellipse.gif'));
            ptImage = ind2rgb(img,map);            
            evalin('base','hlink.Enabled=''off''');
            tt.CData = ptImage;
        end        
    end
end


