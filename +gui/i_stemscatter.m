function [h1, h2] = i_stemscatter(s, c)

c(c < 0) = 0;
x = s(:, 1);
y = s(:, 2);
h1 = stem3(x, y, c, 'marker', 'none', 'color', 'm');
hold on;
h2 = scatter3(x, y, zeros(size(y)), 5, c, 'filled');
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
gui.i_setautumncolor(c, a, true, any(c==0));
