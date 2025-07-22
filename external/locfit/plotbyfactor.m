function plotbyfactor(x, y, f)
%
% very simplistic version of plotbyfactor.
% no symbols, will crash if more than 8 levels.
%

z = unique(f);
cols = ['b', 'r', 'g', 'm', 'c', 'y', 'k', 'w'];

for ix = 1:length(z)
    u = find(f == z(ix));
    plot(x(u), y(u), '.', 'color', cols(ix));
    hold on
end
hold off
