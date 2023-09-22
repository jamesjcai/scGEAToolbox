function [f0] = sc_multiembeddings(s1, s2, ttl1, ttl2)
if nargin < 4, ttl2 = ""; end
if nargin < 3, ttl1 = ""; end
if nargin < 1
    pw1 = fileparts(mfilename('fullpath'));
    load(fullfile(pw1, 'doubleS.mat'), 's1', 's2');
    ttl1 = 'tSNE';
    ttl2 = 'UMAP';
end
f0 = figure('Visible', false);
subplot(1, 2, 1);
if size(s1, 2) == 2
    h1 = scatter(s1(:, 1), s1(:, 2));
    %     xline(0,'r')
    %     yline(0,'r')
else
    h1 = scatter3(s1(:, 1), s1(:, 2), s1(:, 3), 5);
end
if ~isempty(ttl1)
    title(ttl1)
end

subplot(1, 2, 2);
if size(s1, 2) == 2
    h2 = scatter(s2(:, 1), s2(:, 2));
    %     xline(0,'r')
    %     yline(0,'r')
else
    h2 = scatter3(s2(:, 1), s2(:, 2), s2(:, 3), 5);

end
if ~isempty(ttl2)
    title(ttl2)
end
% g3=subplot(2,2,3);
% h3=scatter3(s2(:,1),s2(:,2),s2(:,3),5);

%h=findobj(f0.Children,'type','axes');
%hlink=linkprop(h,{'CameraPosition','CameraUpVector'})
%hlink=linkprop([g1 g2],{'CameraPosition','CameraUpVector'});
evalin('base', 'h=findobj(gcf,''type'',''axes'');');
evalin('base', 'hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');


rotate3d(f0, 'on');
f0.Position(3) = f0.Position(3) * 2;

hBr = brush(f0);
hBr.ActionPostCallback = {@onBrushAction, h1, h2};
movegui(f0, 'center');
set(f0, 'Visible', true);
end

function onBrushAction(~, event, h1, h2)
if isequal(event.Axes.Children, h1)
    h2.BrushData = h1.BrushData;
elseif isequal(event.Axes.Children, h2)
    h1.BrushData = h2.BrushData;
end
% if isprop(h1,'BrushData') && any(h1.BrushData)
%     ptsSelected=logical(h1.BrushData);
%     sum(ptsSelected)
% end
end
