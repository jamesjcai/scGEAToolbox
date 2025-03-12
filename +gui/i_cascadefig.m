function [hFig] = i_cascadefig(sce, g, axx, bxx, ~, methodid)

if nargin < 6, methodid = 5; end
if nargin < 5, k = 1; end

hx=gui.myFigure;
hFig = hx.FigHandle;
[h1] = sc_scattermarker(sce.X, sce.g, sce.s, g, methodid);
view(h1, axx, bxx);
hx.show;

end