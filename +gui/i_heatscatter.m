function [h1] = i_heatscatter(s, c)

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
h1 = scatter3(x, y, z, sz, c, 'filled');
grid on

