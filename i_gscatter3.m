function i_gscatter3(s,c,methodid)

if nargin<3, methodid=1; end
if nargin<2, c='b'; end

x=s(:,1);
y=s(:,2);

if iscell(c), c=grp2idx(c); end
isbinary=logical(numel(unique(c))==2);
        if size(s,2)>=3, z=s(:,3); end

switch methodid
    case 1
        if size(s,2)==2
           scatter(x,y,10,c,'filled');
        elseif size(s,2)>=3
           scatter3(x,y,z,10,c,'filled');
        end
    case 2
        if size(s,2)==2
            gscatter(x,y,c,[],[],10);
        elseif size(s,2)>=3
            gscatter3b(x,y,z,c);
        end
end

if isbinary
    colormap(lines(2));
end

add_3dcamera;

% hFig=gcf;
% tb = uitoolbar(hFig);
% pt = uipushtool(tb,'Separator','off');
% [img,map] = imread(fullfile(fileparts(which(mfilename)),...
%             'private','camera.gif'));
% ptImage = ind2rgb(img,map);
% pt.CData = ptImage;
% pt.Tooltip = 'Select a gene to show expression';
% pt.ClickedCallback = @camera3dmp4;

end
% 
% function camera3dmp4(~,~)
%     OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     fname=tempname;
%     CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],fname,OptionZ)
%     winopen(tempdir);
%     winopen(sprintf('%s.mp4',fname));
% end


function gscatter3b(x,y,z,group,clr,sym,siz,doleg,xnam,ynam,znam)
%GSCATTER3B  3D Scatter plot with grouping variable
%   gscatter3(x,y,z,group,clr,sym,siz,doleg,xnam,ynam,znam)   
%   Designed to work in the exactly same fashion as statistics toolbox's
%   gscatter. Differently from GSCATTER3, this function requires the
%   statistics toolbox installed.
%   
%
%   See also GSCATTER, GSCATTER3
%
%   Copyright 2017 Gustavo Ferraz Trinade.
% Set number of groups 
cgroups = unique(group);
cmap = lines(size(cgroups,1));
% Input variables
if (nargin < 5),  clr = lines(max(size(cgroups))); end
if (nargin < 6) || isempty(sym), sym = '.'; end
if (nargin < 7),  siz = 10;           end
if (nargin < 8),  doleg = 'on';        end
if (nargin < 9),  xnam = inputname(1); end
if (nargin < 10), ynam = inputname(2); end
if (nargin < 11), znam = inputname(3); end
% Get current axes
a = gca;
hold(a,'on')
% call GSCATTER and capture output argument (handles to lines)
h = gscatter(x, y, group,clr,sym,siz,doleg,xnam,ynam);
for i = 1:max(size(cgroups))
        if iscell(cgroups) || ischar(cgroups)
            gi = find(strcmp(group,cgroups(i)));
        else
            gi = find(group == cgroups(i));
        end
    
        set(h(i), 'ZData', z( gi ));
end
zlabel(a,znam);    
view(3)
box on
grid on
    
end