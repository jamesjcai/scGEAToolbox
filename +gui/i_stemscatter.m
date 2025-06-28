function [h1, h2] = i_stemscatter(s, c, ax)

if nargin<3, ax=[]; end

c(c < 0) = 0;
x = s(:, 1);
y = s(:, 2);
if isempty(ax)
    ax = gca; % Get current axes if not provided
end
h1 = stem3(ax, x, y, c, 'marker', 'none', 'color', 'm');
hold(ax,"on");
h2 = scatter3(ax, x, y, zeros(size(y)), 5, c, 'filled');
grid(ax, "on");
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
gui.i_setautumncolor(c, a, true, any(c==0), ax);
