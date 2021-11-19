function [f0]=sc_multigroupings(sce,cx1,cx2,ttl1,ttl2)
if nargin<5, ttl2=""; end
if nargin<4, ttl1=""; end

c1=cx1.c; 
cL1=cx1.cL;

c2=cx2.c;
cL2=cx2.cL;

f0=figure('Visible',false);

ax1=subplot(1,2,1);
h1=gui.i_gscatter3(sce.s,c1,1,1);
if ~isempty(ttl1), title(ax1,ttl1); end
dt = datacursormode(f0);
dt.UpdateFcn = {@i_myupdatefcnx12};

ax2=subplot(1,2,2);
h2=gui.i_gscatter3(sce.s, c2,1,1);
if ~isempty(ttl2), title(ax2,ttl2); end
dt = datacursormode(f0);
dt.UpdateFcn = {@i_myupdatefcnx12};
sgtitle(sce.title);

kc1=numel(unique(c1));
if kc1<=50 && kc1>0
    colormap(ax1,lines(kc1));    
else
    colormap(ax1,'default');
end


kc2=numel(unique(c2));
if kc2<=50 && kc2>0
    colormap(ax2,lines(kc2));    
else
    colormap(ax2,'default');
end

colormap(ax1,lines(kc1));  
evalin('base','h=findobj(gcf,''type'',''axes'');');
evalin('base','hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
rotate3d(f0,'on');
f0.Position(3)=f0.Position(3)*2;

hBr=brush(f0);
hBr.ActionPostCallback = {@onBrushAction,h1,h2};


tb = uitoolbar(f0);

pt = uipushtool(tb, 'Separator', 'off');
[img, map] = imread(fullfile(fileparts(mfilename('fullpath')), ...
                             '..','resources', 'plottypectl-rlocusplot.gif'));  % plotpicker-pie
ptImage = ind2rgb(img, map);
pt.CData = ptImage;
pt.Tooltip = 'Link subplots';
pt.ClickedCallback = @gui.i_linksubplots;
                
gui.add_3dcamera(tb);
movegui(f0,'center');
set(f0,'Visible',true);

function [txt] = i_myupdatefcnx12(Target, event_obj)
    % pos = event_obj.Position;
    if Target.Parent==ax1
        idx = event_obj.DataIndex;
        txt = cL1(c1(idx));
    elseif Target.Parent==ax2
        idx = event_obj.DataIndex;
        txt = cL2(c2(idx));
    end
end

end

function onBrushAction(~,event,h1,h2)
if isequal(event.Axes.Children,h1)
    h2.BrushData=h1.BrushData;
elseif isequal(event.Axes.Children,h2)
    h1.BrushData=h2.BrushData;
end
% if isprop(h1,'BrushData') && any(h1.BrushData)
%     ptsSelected=logical(h1.BrushData);
%     sum(ptsSelected)
% end
end


