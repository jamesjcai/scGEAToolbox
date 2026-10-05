function [h1, h2] = i_stemscatter(s, c, ax)

if nargin<3, ax=[]; end

% Negative values are drawn as downward stems. They used to be clipped
% to zero, so a negative cell-state score looked the same as none.
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
% Grey marks zero only when zero is the bottom of the scale; with negative
% values the lowest colour belongs to the most negative cells.
gui.i_setautumncolor(c, a, true, any(c==0) && min(c) >= 0, ax);
