function [hFig] = i_cascadefig(sce, g, axx, bxx, ~, methodid)

if nargin < 6, methodid = 5; end
if nargin < 5, k = 1; end
hFig = figure('visible', 'off');

[h1] = sc_scattermarker(sce.X, sce.g, sce.s, g, methodid);
view(h1, axx, bxx);
movegui(hFig,'center');
% P = get(f, 'Position');
% set(f, 'Position', [P(1) - 20 * k, P(2) - 20 * k, P(3), P(4)]);
set(hFig, 'visible', 'on');
drawnow;
end