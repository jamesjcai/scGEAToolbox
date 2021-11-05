function [f0]=sc_multigroupings(sce,ttl1,ttl2)
if nargin<3, ttl2=""; end
if nargin<2, ttl1=""; end

if isempty(sce.c_cluster_id) || size(sce.s,1)~=length(sce.c_cluster_id)
    return; 
end
if isempty(sce.c_batch_id)  || size(sce.s,1)~=length(sce.c_batch_id)
    return;
end

f0=figure('Visible',false);

subplot(1,2,1);
h1=gui.i_gscatter3(sce.s, sce.c_cluster_id,1,1);
if ~isempty(ttl1), title(ttl1); end

subplot(1,2,2);
h2=gui.i_gscatter3(sce.s, sce.c_cluster_id,1,1);
if ~isempty(ttl2), title(ttl2); end

evalin('base','h=findobj(gcf,''type'',''axes'');');
evalin('base','hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
rotate3d(f0,'on');
f0.Position(3)=f0.Position(3)*2;

hBr=brush(f0);
hBr.ActionPostCallback = {@onBrushAction,h1,h2};
movegui(f0,'center');
set(f0,'Visible',true);
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
