function sc_multiembeddings(s1,s2)
if nargin<1
    pw1=fileparts(mfilename('fullpath'));
    load(fullfile(pw1,'doubleS.mat'),'s1','s2');
end
f0=figure('Visible',false);
g1=subplot(1,2,1);
h1=scatter3(s1(:,1),s1(:,2),s1(:,3),5);
title('tSNE')
g2=subplot(1,2,2);
h2=scatter3(s2(:,1),s2(:,2),s2(:,3),5);
title('UMAP')
% g3=subplot(2,2,3);
% h3=scatter3(s2(:,1),s2(:,2),s2(:,3),5);

%h=findobj(f0.Children,'type','axes');
%hlink=linkprop(h,{'CameraPosition','CameraUpVector'})
%hlink=linkprop([g1 g2],{'CameraPosition','CameraUpVector'});
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

