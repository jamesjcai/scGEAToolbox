function [h] = i_gscatter3new(s, c, methodid, targetc, hAx)

if nargin < 5, hAx = []; end
if nargin < 4, targetc = 1; end
if nargin < 3, methodid = 1; end
if nargin < 2, c = ones(size(s, 1), 1); end

x = s(:, 1);
y = s(:, 2);

if iscell(c) || isstring(c)
    c = grp2idx(c);
end
kc = numel(unique(c));

if size(s, 2) >= 3, z = s(:, 3); end

switch methodid
    case 1
        % https://www.mathworks.com/matlabcentral/answers/897317-add-legend-to-plot-colored-by-colormap-function
        if size(s, 2) == 2
            if isempty(hAx)                
                h = scatter(x(c==1), y(c==1), 10);
                hold on
                for k=2:max(c)
                    scatter(x(c==k), y(c==k), 10)
                end
            else
                h = scatter(hAx, x(c==1), y(c==1), 10);
                hold on
                for k=2:max(c)
                    scatter(hAx, x(c==k), y(c==k), 10)
                end
            end
        elseif size(s, 2) >= 3
            if isempty(hAx)
                h = scatter3(x(c==1), y(c==1), z(c==1), 10);
                hold on
                for k=2:max(c)
                    scatter3(x(c==k), y(c==k), z(c==k), 10);
                end
            else
                h = scatter3(hAx, x(c==1), y(c==1), z(c==1), 10);
                hold on
                for k=2:max(c)
                    scatter3(hAx, x(c==k), y(c==k), z(c==k), 10);
                end
            end
        end
        %            i=c==1;
        %            h=scatter3(x(i),y(i),z(i),10);
        %            hold on
        %            for k=2:max(c)
        %                i=c==k;
        %                scatter3(x(i),y(i),z(i),10);
        %            end
    case 2
        if size(s, 2) == 2
            h = gscatter(x, y, c, [], [], 5, 'off');
        elseif size(s, 2) >= 3
            h = gscatter3b(x, y, z, c, [], [], 5, 'off');
        end
        box off
    case 3
        h = gui.i_gscatter3(s, c, 1);
        h.MarkerEdgeAlpha = 0;
        hold on
        idx = c == targetc;
        h = gui.i_gscatter3(s(idx, :), c(idx), 1);
        hold off
end

grid on
colormap(pkg.i_mycolorlines(kc));
end
