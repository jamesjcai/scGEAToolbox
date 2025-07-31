function [h] = i_gscatter3(s, c, methodid, targetc, hAx)

if nargin < 5, hAx = []; end
if nargin < 4, targetc = 1; end
if nargin < 3, methodid = 1; end
if nargin < 2, c = ones(size(s, 1), 1); end

x = s(:, 1);
y = s(:, 2);

%if iscell(c) || isstring(c)
    c = findgroups(string(c));
%end
kc = numel(unique(c));

if size(s, 2) >= 3, z = s(:, 3); end

switch methodid
    case 1
        if size(s, 2) == 2
            if isempty(hAx)
                h = scatter(x, y, 10, c);
            else
                h = scatter(hAx, x, y, 10, c);
            end
        elseif size(s, 2) >= 3
            
            if isempty(hAx)
                h = scatter3(x, y, z, 10, c);
            else
                h = scatter3(hAx, x, y, z, 10, c);
            end
        end
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


    if isempty(hAx)
        grid on
        colormap(pkg.i_mycolorlines(kc));
    else
        grid(hAx,'on');
        colormap(hAx,pkg.i_mycolorlines(kc));
    end
end
