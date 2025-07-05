function [h1] = i_heatscatter(s, c, ax)

if nargin<3, ax=[]; end

if size(s, 2) >= 3
    x = s(:, 1);
    y = s(:, 2);
    z = s(:, 3);
    is2d = false;
else
    x = s(:, 1);
    y = s(:, 2);
    z = zeros(size(x));
    is2d = true;
end

% c(c < 0) = 0;
sz = 5;

if is2d
    if isempty(ax)
        h1 = scatter(x, y, sz, c, 'filled');
    else
        h1 = scatter(ax, x, y, sz, c, 'filled');
    end
else
    if isempty(ax)
        h1 = scatter3(x, y, z, sz, c, 'filled');
    else        
        h1 = scatter3(ax, x, y, z, sz, c, 'filled');
    end
end
if isempty(ax)
    grid on
else
    grid(ax, "on");
end

