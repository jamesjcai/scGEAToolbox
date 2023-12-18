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

c(c < 0) = 0;
% x=s(:,1);
% y=s(:,2);
% h1=stem3(x, y, c, 'marker', 'none', 'color', 'm');
% hold on;
% h2=scatter3(x, y, zeros(size(y)), 5, c, 'filled');
% gui.i_setautumncolor(c);
sz = 5;
h1 = scatter3(x, y, z, sz, c, 'filled');
%set(h1,'XTickLabel',[]);
%set(h1,'YTickLabel',[]);
%set(h1,'ZTickLabel',[]);
grid on
% h1.YDataSource='explorer2IDX';
%title(targetg)


% a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
% gui.i_setautumncolor(c, a, false);



%cb=colorbar(h1);
%cb.Label.String =  'Expression Level';
